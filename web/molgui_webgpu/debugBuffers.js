export function _padInt(i, w = 4) {
    const s = String(i | 0);
    return s.length >= w ? s : (' '.repeat(w - s.length) + s);
}

export function dumpVec3Forces(label, arr, n, opts = {}) {
    const { stride = 4, fixed = 6 } = opts;
    if (!arr) throw new Error(`dumpVec3Forces: arr is null for ${label}`);
    const lines = [];
    lines.push(`# ${label} n=${n} stride=${stride} cols=3`);
    for (let i = 0; i < n; i++) {
        const base = i * stride;
        const xs = [];
        xs.push(Number(arr[base + 0]).toFixed(fixed));
        xs.push(Number(arr[base + 1]).toFixed(fixed));
        xs.push(Number(arr[base + 2]).toFixed(fixed));
        lines.push(`${label}[${_padInt(i, 4)}] ${xs.join(' ')}`);
    }
    return lines;
}

export function dumpVec4BufferLines(label, arr, n, opts = {}) {
    const { stride = 4, cols = 3, fixed = 6 } = opts;
    if (!arr) throw new Error(`dumpVec4BufferLines: arr is null for ${label}`);
    const lines = [];
    lines.push(`# ${label} n=${n} stride=${stride} cols=${cols}`);
    for (let i = 0; i < n; i++) {
        const base = i * stride;
        const xs = [];
        for (let k = 0; k < cols; k++) xs.push(Number(arr[base + k]).toFixed(fixed));
        lines.push(`${label}[${_padInt(i, 4)}] ${xs.join(' ')}`);
    }
    return lines;
}

export function dumpAtomParamsLines(label, paramsVec4, n, opts = {}) {
    const { fixed = 6 } = opts;
    if (!paramsVec4) throw new Error(`dumpAtomParamsLines: paramsVec4 is null for ${label}`);
    const lines = [];
    lines.push(`# ${label} n=${n} layout=vec4(radius,0,0,mass)`);
    for (let i = 0; i < n; i++) {
        const base = i * 4;
        const r = Number(paramsVec4[base + 0]).toFixed(fixed);
        const m = Number(paramsVec4[base + 3]).toFixed(fixed);
        lines.push(`${label}[${_padInt(i, 4)}] r=${r} m=${m}`);
    }
    return lines;
}

export function dumpBondIndicesLines(label, bondIdx, nAtoms, nMaxBonded) {
    if (!bondIdx) throw new Error(`dumpBondIndicesLines: bondIdx is null for ${label}`);
    const lines = [];
    lines.push(`# ${label} nAtoms=${nAtoms} nMaxBonded=${nMaxBonded}`);
    for (let i = 0; i < nAtoms; i++) {
        const base = i * nMaxBonded;
        const row = [];
        for (let k = 0; k < nMaxBonded; k++) row.push(String(bondIdx[base + k] | 0));
        lines.push(`${label}[${_padInt(i, 4)}] ${row.join(' ')}`);
    }
    return lines;
}

export function dumpBondLenStiffLines(label, bondLenStiff, nAtoms, nMaxBonded, opts = {}) {
    const { fixed = 6 } = opts;
    if (!bondLenStiff) throw new Error(`dumpBondLenStiffLines: bondLenStiff is null for ${label}`);
    const lines = [];
    lines.push(`# ${label} nAtoms=${nAtoms} nMaxBonded=${nMaxBonded} layout=(L,K) per slot`);
    for (let i = 0; i < nAtoms; i++) {
        const base = i * nMaxBonded;
        const row = [];
        for (let k = 0; k < nMaxBonded; k++) {
            const L = Number(bondLenStiff[(base + k) * 2 + 0]).toFixed(fixed);
            const K = Number(bondLenStiff[(base + k) * 2 + 1]).toFixed(fixed);
            row.push(`${L},${K}`);
        }
        lines.push(`${label}[${_padInt(i, 4)}] ${row.join(' ')}`);
    }
    return lines;
}

export function dumpGhostPackedLines(label, ghostPacked, numGroups, maxGhosts) {
    if (!ghostPacked) throw new Error(`dumpGhostPackedLines: ghostPacked is null for ${label}`);
    const lines = [];
    lines.push(`# ${label} numGroups=${numGroups} maxGhosts=${maxGhosts} layout=[idx... idx][count@maxGhosts]`);
    const stride = maxGhosts + 1;
    for (let g = 0; g < numGroups; g++) {
        const base = g * stride;
        const count = ghostPacked[base + maxGhosts] | 0;
        lines.push(`${label}[grp=${_padInt(g, 3)}] count=${count}`);
        for (let k = 0; k < count; k++) {
            lines.push(`${label}[grp=${_padInt(g, 3)}][${_padInt(k, 3)}] ${ghostPacked[base + k] | 0}`);
        }
    }
    return lines;
}

export function joinSections(sections) {
    const out = [];
    for (const sec of sections) {
        if (!sec) continue;
        for (const l of sec) out.push(l);
    }
    return out.join('\n') + '\n';
}
