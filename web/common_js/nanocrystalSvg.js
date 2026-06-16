/// Orthographic SVG export for nanocrystal debug views (Cartesian-aligned projections).
const Z_COLOR = { 1: '#ffffff', 6: '#404040', 7: '#3050f8', 8: '#ff0d0d', 14: '#daa520' };
const Z_R = { 1: 0.35, 6: 0.77, 7: 0.74, 8: 0.73, 14: 1.11 };

const S2 = 1 / Math.sqrt(2);
const S3 = 1 / Math.sqrt(3);

export const CRYSTAL_VIEWS = {
    '001': { label: '001 (xy)', R: [1, 0, 0, 0, 1, 0, 0, 0, 1] },
    '100': { label: '100 (yz)', R: [0, 1, 0, 0, 0, 1, 1, 0, 0] },
    '010': { label: '010 (xz)', R: [1, 0, 0, 0, 0, 1, 0, 1, 0] },
    '111': { label: '111', R: [S2, -S2, 0, 0, 0, -1, S3, S3, S3] },
};

export function viewMatrixNamed(name = '111') {
    const v = CRYSTAL_VIEWS[name] || CRYSTAL_VIEWS['111'];
    return v.R;
}

export function defaultViewMatrix() { return viewMatrixNamed('111'); }

function mat3MulVec(R, x, y, z) {
    return [R[0] * x + R[1] * y + R[2] * z, R[3] * x + R[4] * y + R[5] * z, R[6] * x + R[7] * y + R[8] * z];
}

function escapeXml(s) {
    return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

function projectCrystal(pos, Z, bonds_ij, R) {
    const n = Z.length;
    const pts = [];
    for (let i = 0; i < n; i++) {
        const [vx, vy] = mat3MulVec(R, pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]);
        pts.push({ i, vx, vy });
    }
    let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    for (const p of pts) {
        xmin = Math.min(xmin, p.vx); xmax = Math.max(xmax, p.vx);
        ymin = Math.min(ymin, p.vy); ymax = Math.max(ymax, p.vy);
    }
    return { pts, xmin, xmax, ymin, ymax, n };
}

function bboxSpan(b) {
    return { spanX: Math.max(b.xmax - b.xmin, 1e-6), spanY: Math.max(b.ymax - b.ymin, 1e-6), cx: 0.5 * (b.xmin + b.xmax), cy: 0.5 * (b.ymin + b.ymax) };
}

function buildLayerSvg(pts, Z, bonds_ij, n, toX, toY, s, bondStroke = '#666') {
    const bondLines = [];
    if (bonds_ij) {
        const nb = bonds_ij.length / 2;
        for (let k = 0; k < nb; k++) {
            const a = bonds_ij[k * 2] | 0, b = bonds_ij[k * 2 + 1] | 0;
            if (a < 0 || b < 0 || a >= n || b >= n) continue;
            bondLines.push(`<line x1="${toX(pts[a].vx).toFixed(2)}" y1="${toY(pts[a].vy).toFixed(2)}" x2="${toX(pts[b].vx).toFixed(2)}" y2="${toY(pts[b].vy).toFixed(2)}" stroke="${bondStroke}" stroke-width="1.2"/>`);
        }
    }
    const sorted = [...pts].sort((a, b) => (Z[a.i] !== Z[b.i] ? Z[a.i] - Z[b.i] : a.vy - b.vy));
    const circles = sorted.map((p) => {
        const rad = (Z_R[Z[p.i]] || 0.8) * s * 0.11;
        const fill = Z_COLOR[Z[p.i]] || '#888';
        const stroke = Z[p.i] === 1 ? '#aaa' : '#222';
        return `<circle cx="${toX(p.vx).toFixed(2)}" cy="${toY(p.vy).toFixed(2)}" r="${Math.max(1.5, rad).toFixed(2)}" fill="${fill}" stroke="${stroke}" stroke-width="0.6"/>`;
    });
    return { bondLines, circles };
}

export function exportCrystalSvg({ pos, Z, bonds_ij = null, width = 520, height = 520, view = '111', R_view = null, title = '', marginFrac = 0.14 }) {
    const R = R_view || viewMatrixNamed(view);
    const viewLabel = (CRYSTAL_VIEWS[view] || CRYSTAL_VIEWS['111']).label;
    const { pts, xmin, xmax, ymin, ymax, n } = projectCrystal(pos, Z, bonds_ij, R);
    const { spanX, spanY, cx, cy } = bboxSpan({ xmin, xmax, ymin, ymax });
    const padX = spanX * marginFrac, padY = spanY * marginFrac;
    const titleH = title ? 22 : 8;
    const marginPx = 28;
    const plotW = width - 2 * marginPx, plotH = height - titleH - marginPx;
    const s = Math.min(plotW / (spanX + 2 * padX), plotH / (spanY + 2 * padY));
    const toX = (x) => width / 2 + (x - cx) * s;
    const toY = (y) => titleH + plotH / 2 - (y - cy) * s;
    const { bondLines, circles } = buildLayerSvg(pts, Z, bonds_ij, n, toX, toY, s);
    const ttl = title ? `<text x="${marginPx}" y="16" font-family="sans-serif" font-size="14" fill="#222">${escapeXml(title)} [${viewLabel}]</text>` : '';
    return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">
  <rect width="100%" height="100%" fill="#f8f8f8"/>
  ${ttl}
  ${bondLines.join('\n  ')}
  ${circles.join('\n  ')}
</svg>`;
}

/// Side-by-side init | relaxed in one SVG (shared view, independent centering per panel).
export function exportCrystalCompareSvg({
    posA, ZA, bondsA, labelA = 'init',
    posB, ZB, bondsB, labelB = 'relaxed',
    width = 1040, height = 520, view = '111', title = '', marginFrac = 0.14, panelGap = 36,
}) {
    const R = viewMatrixNamed(view);
    const viewLabel = (CRYSTAL_VIEWS[view] || CRYSTAL_VIEWS['111']).label;
    const projA = projectCrystal(posA, ZA, bondsA, R);
    const projB = projectCrystal(posB, ZB, bondsB, R);
    const boxA = bboxSpan(projA), boxB = bboxSpan(projB);
    const titleH = title ? 24 : 8;
    const labelH = 18;
    const marginPx = 24;
    const panelW = (width - panelGap - 2 * marginPx) / 2;
    const plotH = height - titleH - labelH - marginPx;
    const padAX = boxA.spanX * marginFrac, padAY = boxA.spanY * marginFrac;
    const padBX = boxB.spanX * marginFrac, padBY = boxB.spanY * marginFrac;
    const sA = Math.min(panelW / (boxA.spanX + 2 * padAX), plotH / (boxA.spanY + 2 * padAY));
    const sB = Math.min(panelW / (boxB.spanX + 2 * padBX), plotH / (boxB.spanY + 2 * padBY));
    const s = Math.min(sA, sB);
    const leftCx = marginPx + panelW / 2;
    const rightCx = marginPx + panelW + panelGap + panelW / 2;
    const midY = titleH + labelH + plotH / 2;
    const toXA = (x) => leftCx + (x - boxA.cx) * s;
    const toYA = (y) => midY - (y - boxA.cy) * s;
    const toXB = (x) => rightCx + (x - boxB.cx) * s;
    const toYB = (y) => midY - (y - boxB.cy) * s;
    const layerA = buildLayerSvg(projA.pts, ZA, bondsA, projA.n, toXA, toYA, s, '#7c6f64');
    const layerB = buildLayerSvg(projB.pts, ZB, bondsB, projB.n, toXB, toYB, s, '#444');
    const ttl = title ? `<text x="${marginPx}" y="16" font-family="sans-serif" font-size="14" fill="#222">${escapeXml(title)} [${viewLabel}]</text>` : '';
    const sepX = marginPx + panelW + panelGap / 2;
    return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">
  <rect width="100%" height="100%" fill="#f8f8f8"/>
  ${ttl}
  <line x1="${sepX.toFixed(1)}" y1="${titleH}" x2="${sepX.toFixed(1)}" y2="${height - 8}" stroke="#ccc" stroke-width="1" stroke-dasharray="4,4"/>
  <text x="${leftCx.toFixed(1)}" y="${titleH + 14}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#8b4513">${escapeXml(labelA)}</text>
  <text x="${rightCx.toFixed(1)}" y="${titleH + 14}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#1d4ed8">${escapeXml(labelB)}</text>
  <g id="panel-init">${layerA.bondLines.join('\n    ')}${layerA.circles.join('\n    ')}</g>
  <g id="panel-relaxed">${layerB.bondLines.join('\n    ')}${layerB.circles.join('\n    ')}</g>
</svg>`;
}

export function exportCrystalCompareSvgViews({ posA, ZA, bondsA, posB, ZB, bondsB, views = ['111', '100', '001'], labelA = 'init', labelB = 'relaxed', title = '', ...rest }) {
    const out = {};
    for (const v of views) {
        out[v] = exportCrystalCompareSvg({ posA, ZA, bondsA, posB, ZB, bondsB, view: v, labelA, labelB, title, ...rest });
    }
    return out;
}

export function exportCrystalSvgViews({ pos, Z, bonds_ij, views = ['111', '100', '001'], title = '', ...rest }) {
    const out = {};
    for (const v of views) out[v] = exportCrystalSvg({ pos, Z, bonds_ij, view: v, title, ...rest });
    return out;
}

/// HTML table row for shape atlas (embeds compare SVG via relative path).
export function atlasTableRow({ id, label, genParams, svgRel, natoms }) {
    const params = escapeXml(JSON.stringify(genParams, null, 0));
    return `<tr><td><b>${escapeXml(label)}</b><br/><code>${escapeXml(id)}</code></td><td><img src="${escapeXml(svgRel)}" width="320" alt="${escapeXml(label)}"/></td><td>${natoms ?? '—'}</td><td><pre style="font-size:10px;max-width:420px;white-space:pre-wrap">${params}</pre></td></tr>`;
}

export function atlasIndexHtml(rows, title = 'Nanocrystal shape atlas') {
    const body = rows.map(atlasTableRow).join('\n');
    return `<!DOCTYPE html><html><head><meta charset="utf-8"/><title>${escapeXml(title)}</title>
<style>body{font-family:sans-serif;margin:20px}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:8px;vertical-align:top}th{background:#eee}</style></head>
<body><h1>${escapeXml(title)}</h1>
<p><a href="viewer.html"><b>Interactive 3D viewer</b></a> (init vs relaxed, drag to rotate)</p>
<table><tr><th>Shape</th><th>Preview [111]</th><th>Atoms</th><th>gen_params</th></tr>
${body}</table></body></html>`;
}
