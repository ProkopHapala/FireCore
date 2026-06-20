/// Numpy .npz/.npy read/write (Node) + nanocrystal crystal geometry NPZ helpers.
import { deflateRawSync, crc32 as zlibCrc32, inflateRawSync } from 'node:zlib';

const Z_TO_SYM = { 1: 'H', 6: 'C', 7: 'N', 8: 'O', 14: 'Si' };

function crc32(buf) {
    if (typeof zlibCrc32 === 'function') return zlibCrc32(buf) >>> 0;
    throw new Error('npzIO: node:zlib.crc32 unavailable; need Node >= 20.12');
}

function descrOf(arr) {
    if (arr instanceof Float64Array) return '<f8';
    if (arr instanceof Float32Array) return '<f4';
    if (arr instanceof Int32Array) return '<i4';
    if (arr instanceof Uint8Array) return '|u1';
    throw new Error(`npzIO: unsupported array type ${arr?.constructor?.name}`);
}

export function encodeNpy(arr, shape) {
    const descr = descrOf(arr);
    const sh = shape && shape.length ? shape : [arr.length];
    const shapeStr = sh.length === 1 ? `(${sh[0]},)` : `(${sh.join(', ')})`;
    let dict = `{'descr': '${descr}', 'fortran_order': False, 'shape': ${shapeStr}, }`;
    let pad = 16 - ((10 + dict.length) % 16);
    if (pad === 16) pad = 0;
    dict += ' '.repeat(pad);
    const hlen = dict.length;
    const header = Buffer.alloc(10 + hlen);
    header.write('\x93NUMPY', 0, 'latin1');
    header.writeUInt8(1, 6);
    header.writeUInt8(0, 7);
    header.writeUInt16LE(hlen, 8);
    header.write(dict, 10, 'latin1');
    const body = Buffer.from(arr.buffer, arr.byteOffset, arr.byteLength);
    return Buffer.concat([header, body]);
}

function zipDeflateEntry(name, raw) {
    const nb = Buffer.from(name, 'utf8');
    const comp = deflateRawSync(raw);
    const crc = crc32(raw);
    const lh = Buffer.alloc(30);
    lh.writeUInt32LE(0x04034b50, 0);
    lh.writeUInt16LE(20, 4);
    lh.writeUInt16LE(0, 6);
    lh.writeUInt16LE(8, 8);
    lh.writeUInt16LE(0, 10);
    lh.writeUInt16LE(0, 12);
    lh.writeUInt32LE(crc, 14);
    lh.writeUInt32LE(comp.length, 18);
    lh.writeUInt32LE(raw.length, 22);
    lh.writeUInt16LE(nb.length, 26);
    lh.writeUInt16LE(0, 28);
    return { local: Buffer.concat([lh, nb, comp]), rawLen: raw.length, compLen: comp.length, crc, nb };
}

function inferShape(key, arr) {
    if (arr.length === 1) return [];
    throw new Error(`npzIO: pass _shape on '${key}' (length=${arr.length})`);
}

export function encodeNpzCompressed(namedArrays) {
    const locals = [];
    const centrals = [];
    let offset = 0;
    for (const [key, arr] of Object.entries(namedArrays)) {
        const name = `${key}.npy`;
        const shape = arr._shape || (arr.length === 1 && !(arr instanceof Float64Array && arr.length > 3) ? [] : null);
        const sh = shape !== null ? shape : inferShape(key, arr);
        const raw = encodeNpy(arr, sh);
        const e = zipDeflateEntry(name, raw);
        e.localOff = offset;
        offset += e.local.length;
        locals.push(e.local);
        const ch = Buffer.alloc(46);
        ch.writeUInt32LE(0x02014b50, 0);
        ch.writeUInt16LE(20, 4);
        ch.writeUInt16LE(20, 6);
        ch.writeUInt16LE(0, 8);
        ch.writeUInt16LE(8, 10);
        ch.writeUInt16LE(0, 12);
        ch.writeUInt16LE(0, 14);
        ch.writeUInt32LE(e.crc, 16);
        ch.writeUInt32LE(e.compLen, 20);
        ch.writeUInt32LE(e.rawLen, 24);
        ch.writeUInt16LE(e.nb.length, 28);
        ch.writeUInt16LE(0, 30);
        ch.writeUInt16LE(0, 32);
        ch.writeUInt16LE(0, 34);
        ch.writeUInt16LE(0, 36);
        ch.writeUInt32LE(0, 38);
        ch.writeUInt32LE(e.localOff, 42);
        centrals.push(Buffer.concat([ch, e.nb]));
    }
    const centralStart = offset;
    const centralBuf = Buffer.concat(centrals);
    const eocd = Buffer.alloc(22);
    eocd.writeUInt32LE(0x06054b50, 0);
    eocd.writeUInt16LE(0, 4);
    eocd.writeUInt16LE(0, 6);
    eocd.writeUInt16LE(Object.keys(namedArrays).length, 8);
    eocd.writeUInt16LE(Object.keys(namedArrays).length, 10);
    eocd.writeUInt32LE(centralBuf.length, 12);
    eocd.writeUInt32LE(centralStart, 16);
    eocd.writeUInt16LE(0, 20);
    return Buffer.concat([...locals, centralBuf, eocd]);
}

export function writeNpzCompressed(path, namedArrays, fs) {
    fs.writeFileSync(path, encodeNpzCompressed(namedArrays));
}

function parseNpy(buf) {
    if (buf.length < 10 || buf.toString('latin1', 0, 6) !== '\x93NUMPY') throw new Error('npzIO: not npy');
    const major = buf.readUInt8(6);
    const hlen = major === 1 ? buf.readUInt16LE(8) : buf.readUInt32LE(8);
    const hoff = major === 1 ? 10 : 12;
    const header = buf.toString('latin1', hoff, hoff + hlen);
    const m = header.match(/'descr':\s*'([^']+)'/);
    const sm = header.match(/'shape':\s*\(([^)]*)\)/);
    if (!m || !sm) throw new Error(`npzIO: bad header: ${header.slice(0, 120)}`);
    const descr = m[1];
    const shapeParts = sm[1].trim().split(',').map(s => s.trim()).filter(Boolean);
    const shape = shapeParts.length ? shapeParts.map(x => parseInt(x, 10)) : [];
    const dataOff = hoff + hlen;
    const raw = buf.subarray(dataOff);
    if (descr === '<f8') return { data: new Float64Array(raw.buffer, raw.byteOffset, raw.byteLength / 8), shape, descr };
    if (descr === '<f4') return { data: new Float32Array(raw.buffer, raw.byteOffset, raw.byteLength / 4), shape, descr };
    if (descr === '<i4') return { data: new Int32Array(raw.buffer, raw.byteOffset, raw.byteLength / 4), shape, descr };
    if (descr === '|u1') return { data: new Uint8Array(raw.buffer, raw.byteOffset, raw.byteLength), shape, descr };
    throw new Error(`npzIO: unsupported descr ${descr}`);
}

function readZipEntries(buf) {
    const out = new Map();
    let off = 0;
    while (off + 30 <= buf.length) {
        const sig = buf.readUInt32LE(off);
        if (sig === 0x06054b50 || sig === 0x02014b50) break;
        if (sig !== 0x04034b50) throw new Error(`npzIO: bad zip sig at ${off}`);
        const compM = buf.readUInt16LE(off + 8);
        const compSz = buf.readUInt32LE(off + 18);
        const rawSz = buf.readUInt32LE(off + 22);
        const nameLen = buf.readUInt16LE(off + 26);
        const extraLen = buf.readUInt16LE(off + 28);
        const name = buf.toString('utf8', off + 30, off + 30 + nameLen);
        const dataStart = off + 30 + nameLen + extraLen;
        const comp = buf.subarray(dataStart, dataStart + compSz);
        let raw;
        if (compM === 0) raw = comp;
        else if (compM === 8) raw = inflateRawSync(comp);
        else throw new Error(`npzIO: unsupported compression ${compM} in ${name}`);
        if (raw.length !== rawSz) throw new Error(`npzIO: size mismatch ${name}`);
        out.set(name.replace(/\.npy$/, ''), parseNpy(raw));
        off = dataStart + compSz;
    }
    return out;
}

export function readNpzFile(fs, filePath) {
    const buf = fs.readFileSync(filePath);
    const entries = readZipEntries(buf);
    const arrays = {};
    for (const [k, v] of entries) arrays[k] = v.data;
    return { arrays, entries };
}

export function readNpzUtf8(entries, key) {
    const e = entries.get(key);
    if (!e) return null;
    const u8 = e.data instanceof Uint8Array ? e.data : new Uint8Array(e.data.buffer, e.data.byteOffset, e.data.byteLength);
    return new TextDecoder('utf-8').decode(u8);
}

export function molToCrystalArrays(mol) {
    const n = mol.atoms.length;
    const pos = new Float64Array(n * 3);
    const Z = new Int32Array(n);
    const bonds = [];
    for (let i = 0; i < n; i++) {
        const a = mol.atoms[i];
        pos[i * 3] = +a.pos.x;
        pos[i * 3 + 1] = +a.pos.y;
        pos[i * 3 + 2] = +a.pos.z;
        Z[i] = a.Z | 0;
    }
    for (const b of mol.bonds || []) {
        b.ensureIndices(mol);
        bonds.push([b.a | 0, b.b | 0]);
    }
    const bonds_ij = new Int32Array(bonds.length * 2);
    for (let k = 0; k < bonds.length; k++) {
        bonds_ij[k * 2] = bonds[k][0];
        bonds_ij[k * 2 + 1] = bonds[k][1];
    }
    return { pos, Z, bonds_ij, natoms: n, nbonds: bonds.length };
}

export function writeCrystalNpz(fs, outPath, { pos, Z, bonds_ij, gen_params = '', timing_ms = 0 }) {
    const n = Z.length;
    pos._shape = [n, 3];
    Z._shape = [n];
    if (bonds_ij && bonds_ij.length) bonds_ij._shape = [bonds_ij.length / 2, 2];
    const genBytes = new Uint8Array(new TextEncoder().encode(String(gen_params)));
    genBytes._shape = [genBytes.length];
    const arrays = { pos, Z, natoms: new Int32Array([n]), gen_params: genBytes };
    if (bonds_ij && bonds_ij.length) arrays.bonds_ij = bonds_ij;
    if (timing_ms > 0) arrays.timing_ms = new Float64Array([timing_ms]);
    writeNpzCompressed(outPath, arrays, fs);
}

export function readCrystalNpz(fs, filePath) {
    const { arrays, entries } = readNpzFile(fs, filePath);
    const n = arrays.natoms ? arrays.natoms[0] : (arrays.Z ? arrays.Z.length : 0);
    const gen_params = readNpzUtf8(entries, 'gen_params') || '';
    return { pos: arrays.pos, Z: arrays.Z, bonds_ij: arrays.bonds_ij || null, natoms: n, gen_params, entries };
}

export function readXyzPositions(fs, xyzPath) {
    const text = fs.readFileSync(xyzPath, 'utf8');
    const lines = text.trim().split('\n');
    const n = parseInt(lines[0], 10);
    const symMap = { H: 1, C: 6, Si: 14 };
    const pos = new Float64Array(n * 3);
    const Z = new Int32Array(n);
    for (let i = 0; i < n; i++) {
        const p = lines[i + 2].trim().split(/\s+/);
        Z[i] = symMap[p[0]] || 6;
        pos[i * 3] = +p[1];
        pos[i * 3 + 1] = +p[2];
        pos[i * 3 + 2] = +p[3];
    }
    return { pos, Z, natoms: n };
}

export function crystalToXYZ(pos, Z) {
    const n = Z.length;
    const lines = [`${n}`, 'nanocrystal pipeline'];
    for (let i = 0; i < n; i++) {
        const sym = Z_TO_SYM[Z[i]] || 'X';
        lines.push(`${sym} ${pos[i * 3].toFixed(6)} ${pos[i * 3 + 1].toFixed(6)} ${pos[i * 3 + 2].toFixed(6)}`);
    }
    return lines.join('\n') + '\n';
}

/// Compact JSON for nanocrystalViewer.html (pos centered at centroid).
export function crystalToJson({ id = '', label = '', stage = '', pos, Z, bonds_ij = null, extra = null }) {
    const n = Z.length;
    let cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < n; i++) { cx += pos[i * 3]; cy += pos[i * 3 + 1]; cz += pos[i * 3 + 2]; }
    cx /= n; cy /= n; cz /= n;
    const posC = new Array(n * 3);
    for (let i = 0; i < n; i++) {
        posC[i * 3] = pos[i * 3] - cx;
        posC[i * 3 + 1] = pos[i * 3 + 1] - cy;
        posC[i * 3 + 2] = pos[i * 3 + 2] - cz;
    }
    const out = {
        id, label, stage, natoms: n,
        Z: Array.from(Z),
        pos: posC,
        bonds_ij: bonds_ij ? Array.from(bonds_ij) : [],
    };

    if (extra) {
        for (const [k, v] of Object.entries(extra)) {
            if (out[k] !== undefined) throw new Error(`crystalToJson: extra key collides with base key '${k}'`);
            out[k] = v;
        }
        if (out.group_bbox_min && out.group_bbox_max) {
            const mn = out.group_bbox_min;
            const mx = out.group_bbox_max;
            if (!Array.isArray(mn) || !Array.isArray(mx) || (mn.length !== mx.length) || (mn.length % 3 !== 0)) throw new Error('crystalToJson: group_bbox_min/max must be flat arrays of length 3*n');
            for (let i = 0; i < mn.length; i += 3) {
                mn[i + 0] -= cx; mn[i + 1] -= cy; mn[i + 2] -= cz;
                mx[i + 0] -= cx; mx[i + 1] -= cy; mx[i + 2] -= cz;
            }
        }
    }

    return out;
}

export function writeCrystalJson(fs, filePath, data) {
    fs.writeFileSync(filePath, JSON.stringify(data));
}
