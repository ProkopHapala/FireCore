
import { Vec3 } from '../../common_js/Vec3.js';
import { Mat3 } from '../../common_js/Mat3.js';
import * as CrystalUtils from './CrystalUtils.js';
import { EditableMolecule } from './EditableMolecule.js';
import { BuildersGUI } from './BuildersGUI.js';

/**
 * ScriptRunner: High-level command dispatcher for MolGUI.
 * No eval/Function used. Whitelist-based command execution from JSON/Text.
 */
export class ScriptRunner {
    constructor(app) {
        this.app = app;
        this.activeSystemName = 'molecule';
        this.commands = {
            'clear':           this.clear.bind(this),
            'new_system':      this.newSystem.bind(this),
            'use_system':      this.useSystem.bind(this),
            'merge_systems':   this.mergeSystems.bind(this),
            'load_molecule':   this.loadMolecule.bind(this),
            'make_substrate':  this.buildSubstrate.bind(this),
            'translate':       this.translate.bind(this),
            'rotate':          this.rotate.bind(this),
            'select':          this.select.bind(this),
            'select_query':    this.selectQuery.bind(this),
            'select_all':      this.selectAll.bind(this),
            'clear_selection': this.clearSelection.bind(this),
            'replicate':       this.replicate.bind(this),
            'bake':            this.bake.bind(this),
            'build_substrate': this.buildSubstrate.bind(this),
            'build_polymer':   this.buildPolymer.bind(this),
            'replication':     this.replication.bind(this)
        };
    }

    get system() {
        const s = this.app.systems[this.activeSystemName];
        return s || this.app.systems.molecule || this.app.system;
    }

    async run(script) {
        if (typeof script === 'string' && !script.trim().startsWith('[') && !script.trim().startsWith('{')) {
            // Execute JavaScript script with a safe, whitelisted API
            await this.runJS(script);
            this.finalize();
            return;
        }
        // Fallback for array/object command lists (legacy JSON)
        let cmds = [];
        try {
            cmds = (typeof script === 'string') ? JSON.parse(script) : script;
            if (!Array.isArray(cmds)) cmds = [cmds];
        } catch (e) {
            console.error("ScriptRunner.run(): failed to parse JSON script:", e);
            if (window.logger) window.logger.error("Script JSON Error: " + e.message);
            return;
        }
        for (const cmdObj of cmds) {
            const handler = this.commands[cmdObj.cmd];
            if (handler) {
                try {
                    await handler(cmdObj);
                } catch (e) {
                    console.error(`ScriptRunner: command "${cmdObj.cmd}" failed:`, e);
                    if (window.logger) window.logger.error(`Command "${cmdObj.cmd}" error: ${e.message}`);
                }
            } else {
                console.warn(`ScriptRunner: unknown command "${cmdObj.cmd}"`);
                if (window.logger) window.logger.warn(`Unknown command: ${cmdObj.cmd}`);
            }
        }
        this.finalize();
    }

    async runJS(text) {
        const queue = [];
        const api = {};
        for (const name of Object.keys(this.commands)) {
            api[name] = (opts = {}) => { queue.push({ name, opts }); };
        }
        api.run = (cmd, opts = {}) => { queue.push({ name: cmd, opts }); };

        // Handle-based system API
        const createSystemHandle = (name) => {
            const sys = this.app.systems[name];
            if (!sys) throw new Error(`createSystemHandle: system '${name}' not found`);
            const rend = this.app.renderers[name];
            const handle = {
                id: name,
                system: sys,
                renderer: rend,
                clear: () => { queue.push({ name: 'clear', opts: { system: name } }); },
                load: (opts) => { 
                    const o = (typeof opts === 'string') ? { path: opts } : opts;
                    queue.push({ name: 'load_molecule', opts: { ...o, system: name } }); 
                },
                build_substrate: (preset, opts) => {
                    let o = (typeof preset === 'string') ? { preset, ...opts } : preset;
                    if (o.size) { o.nx = o.size[0]; o.ny = o.size[1]; o.nz = o.size[2]; }
                    queue.push({ name: 'build_substrate', opts: { ...o, system: name } });
                },
                move: (vec) => { queue.push({ name: 'translate', opts: { vec, system: name } }); },
                translate: (vec) => { queue.push({ name: 'translate', opts: { vec, system: name } }); },
                rotate: (axis, deg, ctr) => { queue.push({ name: 'rotate', opts: { axis, deg, center: ctr, system: name } }); },
                roate: (axis, deg, ctr) => { queue.push({ name: 'rotate', opts: { axis, deg, center: ctr, system: name } }); },
                replicate: (nrep) => { queue.push({ name: 'replicate', opts: { n: nrep, system: name } }); },
                replication: (opts) => { queue.push({ name: 'replication', opts: { ...opts, name: name, system: name } }); }
            };
            return handle;
        };

        api.get_system = createSystemHandle;
        const molecule  = createSystemHandle('molecule');
        const substrate = createSystemHandle('substrate');
        const mol = molecule;

        const decls = Object.keys(this.commands).map(n => `const ${n} = api.${n};`).join('\n');

        const fn = new Function('api', 'app', 'CrystalUtils', 'EditableMolecule', 'Vec3', 'Mat3', 'console', 'logger', 'molecule', 'substrate', 'mol', `
            "use strict";
            ${decls}
            const get_system = api.get_system;
            ${text}
        `);
        try {
            fn(api, this.app, CrystalUtils, EditableMolecule, Vec3, Mat3, console, window.logger || console, molecule, substrate, mol);
            for (const { name, opts } of queue) {
                const cmd = this.commands[name];
                if (!cmd) throw new Error(`Unknown command '${name}'`);
                await cmd(opts || {});
            }
        } catch (e) {
            console.error('ScriptRunner.runJS error:', e);
            if (window.logger) window.logger.error(String(e));
        }
    }

    finalize() {
        if (!this.app) return;
        Object.values(this.app.renderers).forEach(r => r.update());
        if (this.app.editor) this.app.editor.updateGizmo();
        if (this.app.requestRender) this.app.requestRender();
    }

    // --- Command Handlers ---

    clear(args = {}) {
        const prev = this.activeSystemName;
        if (args.system) this.useSystem(args.system);
        this.system.clear();
        if (args.system) this.useSystem(prev);
    }

    newSystem(args) {
        const name = (typeof args === 'string') ? args : (args.name || 'new');
        if (!this.app.systems[name]) {
            this.app.systems[name] = new EditableMolecule();
            const packed = new PackedMolecule();
            this.app.renderers[name] = new MoleculeRenderer(this.app.scene, packed, this.app.shaders, this.app.mmParams, this.app.systems[name]);
            if (this.app.renderers[name].toggleAxes) this.app.renderers[name].toggleAxes(true);
        }
        this.activeSystemName = name;
        return this.app.systems[name];
    }

    useSystem(args) {
        const name = (typeof args === 'string') ? args : (args.name || 'molecule');
        if (!this.app.systems[name]) {
            return this.newSystem(name);
        }
        this.activeSystemName = name;
    }

    mergeSystems(args) {
        const { source, target = 'molecule', deleteSource = true } = args;
        if (!source) throw new Error('merge_systems: source name required');
        const src = this.app.systems[source];
        const dst = this.app.systems[target];
        if (!src || !dst) throw new Error(`merge_systems: system missing src=${!!src} dst=${!!dst}`);

        const parsed = src.exportAsParsed();
        dst.appendParsedSystem(parsed);
        
        if (deleteSource && source !== 'main' && source !== 'molecule' && source !== 'substrate') {
            delete this.app.systems[source];
            delete this.app.renderers[source];
            if (this.activeSystemName === source) this.activeSystemName = target;
        }
        dst.dirtyExport = true;
    }

    async loadMolecule(args) {
        const { path, fmt, system } = args;
        if (!path) throw new Error("load_molecule: path required");
        
        const prev = this.activeSystemName;
        if (system) this.useSystem(system);

        const realFmt = fmt || path.split('.').pop().toLowerCase();
        
        const response = await fetch(path);
        if (!response.ok) throw new Error(`load_molecule: failed to fetch ${path}`);
        const text = await response.text();
        
        const parsed = (realFmt === 'mol' || realFmt === 'mol2')
            ? EditableMolecule.parseMol2(text)
            : EditableMolecule.parseXYZ(text);

        const mm = this.system;
        const n0 = mm.atoms.length;
        const ids = mm.appendParsedSystem(parsed);
        // Auto-bond XYZ (mol/mol2 already have bonds)
        if (realFmt === 'xyz') {
            const mp = this.app.mmParams || null;
            mm.recalculateBonds(mp);
        }
        mm.clearSelection();
        ids.forEach(id => mm.selection.add(id));
        mm.dirtyExport = true;

        if (system) this.useSystem(prev);
    }

    makeSubstrate(args) {
        this.buildSubstrate(args);
    }

    translate(args) {
        const { vec, sel = 'current', system } = args;
        if (!vec || vec.length < 3) throw new Error("translate: vec [x,y,z] required");
        
        const prev = this.activeSystemName;
        if (system) this.useSystem(system);

        const mm = this.system;
        const d = (vec instanceof Vec3) ? vec : new Vec3(vec[0], vec[1], vec[2]);
        const ids = (sel === 'all') ? mm.atoms.map(a => a.id) : Array.from(mm.selection);
        mm.translateAtoms(ids, d);

        if (system) this.useSystem(prev);
    }

    rotate(args) {
        const { axis, deg, sel = 'current', center, system } = args;
        if (!axis || axis.length < 3 || deg === undefined) throw new Error("rotate: axis and deg required");
        
        const prev = this.activeSystemName;
        if (system) this.useSystem(system);

        const mm = this.system;
        const ids = (sel === 'all') ? mm.atoms.map(a => a.id) : Array.from(mm.selection);
        const vAxis = (axis instanceof Vec3) ? axis : new Vec3(axis[0], axis[1], axis[2]);
        const ctr = (center && center.length >= 3) ? new Vec3(center[0], center[1], center[2]) : null;
        mm.rotateAtoms(ids, vAxis, deg, ctr);

        if (system) this.useSystem(prev);
    }

    select(args) {
        const { ids, mode = 'replace' } = args;
        if (!ids) throw new Error("select: ids required");
        const mm = this.system;
        if (mode === 'replace') mm.selection.clear();
        ids.forEach(id => {
            if (mm.id2atom.has(id)) mm.selection.add(id);
        });
        mm.dirtyExport = true;
    }

    selectQuery(args) {
        const { query, mode = 'replace', print = false } = args;
        if (!query || !window.app || !window.app.mmParams) throw new Error("select_query: requires query and mmParams");
        const cq = EditableMolecule.compileSelectQuery(query, window.app.mmParams);
        this.system.applySelectQuery(cq, { mode, bPrint: !!print });
        this.system.dirtyExport = true;
    }

    selectAll() {
        this.system.selectAll();
    }

    clearSelection() {
        this.system.clearSelection();
    }

    replicate(args) {
        const system = args.system || 'default';
        this.replication({ name: system, system, ...args });
    }

    bake() {
        this.app.bakeReplicas();
    }

    replication(args) {
        const { name = 'default', n, lattice, show = true, showBox = false, system } = args;
        const lat = this.app.getLattice(name);
        if (n) { lat.nrep.x = n[0]; lat.nrep.y = n[1]; lat.nrep.z = n[2]; }
        
        // Sync lattice vectors from system if requested
        const targetSys = system ? this.app.systems[system] : this.app.systems[name];
        if (targetSys && targetSys.lvec && !lattice) {
            lat.lvec[0].setV(targetSys.lvec[0]);
            lat.lvec[1].setV(targetSys.lvec[1]);
            lat.lvec[2].setV(targetSys.lvec[2]);
        }

        if (lattice) {
            lat.lvec[0].set(lattice[0][0], lattice[0][1], lattice[0][2]);
            lat.lvec[1].set(lattice[1][0], lattice[1][1], lattice[1][2]);
            lat.lvec[2].set(lattice[2][0], lattice[2][1], lattice[2][2]);
        }
        lat.show = show;
        lat.showBox = showBox;
        this.app.updateReplicas(name);
    }

    // --- Builder bindings ---

    buildSubstrate(args) {
        const { a = 2.82065, nx = 13, ny = 12, nz = 3, replication = true, system, preset, step_edge } = args || {};
        
        const prev = this.activeSystemName;
        if (system) this.useSystem(system);

        // TODO: Handle 'preset' and 'step_edge' properly if they differ from default genNaClStep
        const mol = CrystalUtils.genNaClStep({ a, nx, ny, nz });
        const parsed = mol.exportAsParsed();
        const ids = this.system.appendParsedSystem(parsed, { pos: new Vec3(0, 0, 0) });
        this.system.clearSelection();
        ids.forEach(id => this.system.selection.add(id));
        this.system.dirtyExport = true;

        if (replication) {
            this.replication({
                name: system || 'substrate',
                system: system || 'substrate',
                n: [0, 0, 0], // Default to 1x1x1 (no extra replicas)
                lattice: [
                    [a * nx, 0, 0],
                    [0, a * ny, 0],
                    [0, 0, a * nz]
                ],
                show: true
            });
        }

        if (system) this.useSystem(prev);
    }

    buildPolymer(args) {
        const { sequence } = args || {};
        if (!sequence) throw new Error('build_polymer: sequence required');
        if (!this.app.gui || !this.app.gui.monomerLib) throw new Error('build_polymer: monomerLib missing (load via GUI first)');
        if (!this.app.gui.monomerGeoms || this.app.gui.monomerGeoms.size === 0) throw new Error('build_polymer: monomer geometries missing (load via GUI)');

        const tokens = this.app.gui.parseSequenceTokens(sequence);
        const monomers = {};
        for (const m of this.app.gui.monomerLib.monomers) {
            const parsed = this.app.gui.monomerGeoms.get(m.file);
            if (!parsed) throw new Error(`build_polymer: missing geometry for monomer '${m.id}' file='${m.file}'`);
            if (!m.anchors || !m.anchors.head || !m.anchors.tail) throw new Error(`build_polymer: monomer '${m.id}' missing anchors.head/tail`);
            if (m.anchors.head.type !== 'index' || m.anchors.tail.type !== 'index') throw new Error(`build_polymer: monomer '${m.id}' anchors must be type='index'`);
            monomers[m.id] = { parsed, anchors: [m.anchors.head.value, m.anchors.tail.value] };
            if (m.aliases) for (const a of m.aliases) monomers[a] = monomers[m.id];
        }
        const poly = BuildersGUI.assemblePolymerFromTokens ? BuildersGUI.assemblePolymerFromTokens(tokens, monomers, { _0: 1 }) : null;
        if (!poly || !poly.atoms || !poly.bonds) throw new Error('build_polymer: assemblePolymerFromTokens unavailable or failed');

        const idsNew = new Map();
        if (args.mode === 'replace') this.system.clear();
        for (const a of poly.atoms) {
            const id = this.system.addAtom(a.pos.x, a.pos.y, a.pos.z, a.Z);
            idsNew.set(a.id, id);
        }
        for (const b of poly.bonds) {
            const aId = idsNew.get(b.aId);
            const cId = idsNew.get(b.bId);
            if (aId === undefined || cId === undefined) throw new Error('build_polymer: bond endpoint missing after copy');
            this.system.addBond(aId, cId);
        }
        this.system.clearSelection();
        idsNew.forEach((id) => this.system.selection.add(id));
        this.system.dirtyExport = true;
    }
}
