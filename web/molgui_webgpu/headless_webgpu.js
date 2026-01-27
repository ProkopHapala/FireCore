import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const _originalFetch = globalThis.fetch;

async function _fileFetch(url) {
    const u = String(url);
    if (u.startsWith('http://') || u.startsWith('https://')) {
        throw new Error(`_fileFetch: HTTP(S) not supported in headless mode url='${u}'`);
    }
    const p = path.resolve(process.cwd(), u);
    const text = await fs.readFile(p, 'utf8');
    return { text: async () => text, ok: true, status: 200 };
}

export async function bootstrap({ dawnModule = 'webgpu' } = {}) {
    if (globalThis.navigator && globalThis.navigator.gpu) return globalThis.navigator.gpu;

    let mod = null;
    try {
        mod = await import(dawnModule);
    } catch (e) {
        const msg = [
            'Node WebGPU backend missing.',
            `Tried to import: ${dawnModule}`,
            '',
            'Install a Dawn WebGPU binding, e.g.:',
            `  npm install --save-dev ${dawnModule} @webgpu/types`,
            '',
            `Original error: ${e?.message || e}`,
        ].join('\n');
        throw new Error(msg);
    }

    const create = mod.create || mod.default?.create;
    const globals = mod.globals || mod.default?.globals;

    if (!create || typeof create !== 'function') {
        throw new Error(`bootstrap: imported '${dawnModule}' but did not find a 'create' function`);
    }

    const gpu = create([]);

    if (!globals) {
        throw new Error(`bootstrap: imported '${dawnModule}' but did not find 'globals' to inject WebGPU constants`);
    }
    Object.assign(globalThis, globals);

    if (!globalThis.navigator) globalThis.navigator = {};
    globalThis.navigator.gpu = gpu;

    const _fetchPolyfill = async (url, options) => {
        const u = String(url);
        if (u.startsWith('http://') || u.startsWith('https://')) {
            throw new Error(`bootstrap.fetch: HTTP(S) not supported in headless mode url='${u}'`);
        }
        const p = path.resolve(path.dirname(fileURLToPath(import.meta.url)), u);
        const text = await fs.readFile(p, 'utf8');
        return { text: async () => text, ok: true, status: 200 };
    };

    Object.defineProperty(globalThis, 'fetch', {
        value: _fetchPolyfill,
        writable: true,
        configurable: true,
        enumerable: true
    });

    return gpu;
}
