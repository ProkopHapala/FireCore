
import { Vec3 } from '../../common_js/Vec3.js';
import { Mat3 } from '../../common_js/Mat3.js';
import * as CrystalUtils from './CrystalUtils.js';
import { EditableMolecule } from './EditableMolecule.js';
import { BuildersGUI } from './BuildersGUI.js';
import { selectBridgeCandidates } from './MoleculeSelection.js';
import { collapseBridgeAt, collapseBridgeRandom, collapseAllBridges, insertBridge, insertBridgeRandom } from './MoleculeUtils.js';

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
            'replication':     this.replication.bind(this),
            'collapse_bridge': this.collapseBridge.bind(this),
            'collapse_bridge_at': this.collapseBridgeAt.bind(this),
            'pick_random_selection': this.pickRandomSelection.bind(this),
            'select_bridge_candidates': this.selectBridgeCandidates.bind(this),
            'collapse_all_bridges': this.collapseAllBridges.bind(this),
            'insert_bridge': this.insertBridge.bind(this),
            'insert_bridge_random': this.insertBridgeRandom.bind(this),
            'addLvec': this.addLvec.bind(this),
            'setViewReplicas': this.setViewReplicas.bind(this)
        };
        this.lastPickedSelection = null;
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
                replication: (opts) => { queue.push({ name: 'replication', opts: { ...opts, name: name, system: name } }); },
                addLvec: (...vecs) => {
                    let lvecParam = vecs;
                    if (vecs.length === 1) {
                        const arg0 = vecs[0];
                        if (arg0 && Array.isArray(arg0) && arg0.length === 3 && !Array.isArray(arg0[0]) && typeof arg0[0] === 'number') {
                            // Single flat numeric array -> treat as full lattice (convert to 3 Vec3 later)
                            lvecParam = [[arg0[0], arg0[1], arg0[2]], [0, 0, 0], [0, 0, 0]];
                        } else if (arg0 && arg0.lvec) {
                            lvecParam = arg0.lvec;
                        } else if (Array.isArray(arg0) && arg0.length === 3) {
                            lvecParam = arg0;
                        }
                    }
                    queue.push({ name: 'addLvec', opts: { lvec: lvecParam, system: name } });
                },
                setViewReplicas: (opts) => { queue.push({ name: 'setViewReplicas', opts: { ...opts, system: name } }); }
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
        let ids = (sel === 'all') ? mm.atoms.map(a => a.id) : Array.from(mm.selection);
        if (ids.length === 0) ids = mm.atoms.map(a => a.id);

        const logFn = (window.logger && typeof window.logger.info === 'function') ? window.logger.info.bind(window.logger) : console.log;
        logFn(`translate: sys=${this.activeSystemName} ids=${ids.length} vec=${d}`);

        mm.translateAtoms(ids, d);
        
        mm.dirtyGeom = true;
        mm.dirtyExport = true;
        const rend = this.app.renderers[this.activeSystemName];
        if (rend) rend.update();
        if (this.app && this.app.requestRender) this.app.requestRender();
        if (system) this.useSystem(prev);
    }

    rotate(args) {
        const { axis, deg, sel = 'current', center, system } = args;
        if (!axis || deg === undefined) throw new Error("rotate: axis and deg required");
        const prev = this.activeSystemName;
        if (system) this.useSystem(system);
        const mm = this.system;
        
        const vAxis = (axis instanceof Vec3) ? axis : new Vec3(axis[0], axis[1], axis[2]);
        let ids = (sel === 'all') ? mm.atoms.map(a => a.id) : Array.from(mm.selection);
        if (ids.length === 0) ids = mm.atoms.map(a => a.id);

        let vCenter;
        if (center) {
            vCenter = (center instanceof Vec3) ? center : new Vec3(center[0], center[1], center[2]);
        } else {
            const points = ids.map(id => {
                const ia = mm.id2atom.get(id);
                return (ia !== undefined) ? mm.atoms[ia].pos : null;
            }).filter(p => p !== null);
            vCenter = Vec3.getCOG(points);
        }

        const logFn = (window.logger && typeof window.logger.info === 'function') ? window.logger.info.bind(window.logger) : console.log;
        logFn(`rotate: sys=${this.activeSystemName} ids=${ids.length} axis=${vAxis} deg=${deg} center=${vCenter}`);

        mm.rotateAtoms(ids, vAxis, deg, vCenter);

        mm.dirtyGeom = true;
        mm.dirtyExport = true;
        const rend = this.app.renderers[this.activeSystemName];
        if (rend) rend.update();
        if (this.app && this.app.requestRender) this.app.requestRender();
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
        // Ensure lattice object exists even if caller runs before init wiring
        if (!this.app.lattices) this.app.lattices = new Map();
        if (!this.app.lattices.has(name)) {
            this.app.lattices.set(name, {
                lvec: [new Vec3(10, 0, 0), new Vec3(0, 10, 0), new Vec3(0, 0, 10)],
                nrep: { x: 0, y: 0, z: 0 },
                show: false,
                showBox: false,
                filter: null
            });
        }
        let lat = this.app.lattices.get(name);
        if (!lat && this.app.getLattice) lat = this.app.getLattice(name);
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
        window.app.updateReplicas(name);
    }

    addLvec(args) {
        const { lvec, system } = args;
        const sysName = system || this.activeSystemName;
        const sys = this.app.systems[sysName];
        if (!sys) return;
        if (lvec && lvec.length === 3) {
            const parsed = [
                this._vecFromInput(lvec[0]),
                this._vecFromInput(lvec[1]),
                this._vecFromInput(lvec[2])
            ];
            sys.lvec = [parsed[0].clone(), parsed[1].clone(), parsed[2].clone()];
            const lat = this.app.getLattice(sysName);
            lat.lvec = [sys.lvec[0].clone(), sys.lvec[1].clone(), sys.lvec[2].clone()];
            const rend = this.app.renderers[sysName];
            if (rend) {
                rend.setReplicas({ lvec: sys.lvec });
                rend.update();
            }
            this.app.requestRender();
        }
    }
    
    setViewReplicas(args) {
        console.log('setViewReplicas', args);
        const { show, showBox, nrep, lvec, system } = args;
        const sysName = system || this.activeSystemName;
        const lat = this.app.getLattice(sysName);
        if (show !== undefined) lat.show = show;
        if (showBox !== undefined) lat.showBox = showBox;
        if (nrep !== undefined) {
            if (Array.isArray(nrep)) { lat.nrep = { x: nrep[0], y: nrep[1], z: nrep[2] }; } 
            else { lat.nrep = { ...lat.nrep, ...nrep }; }
        }
        if (lvec) {
            lat.lvec = [
                this._vecFromInput(lvec[0]),
                this._vecFromInput(lvec[1]),
                this._vecFromInput(lvec[2])
            ];
        }
        console.log('setViewReplicas lvec', lat.lvec[0], lat.lvec[1], lat.lvec[2]);
        this.app.updateReplicas(sysName);
    }

    _vecFromInput(v) {
        if (!v) return new Vec3();
        if (v instanceof Vec3) return v.clone();
        if (Array.isArray(v)) return new Vec3(v[0] || 0, v[1] || 0, v[2] || 0);
        if (typeof v === 'object' && 'x' in v && 'y' in v && 'z' in v) return new Vec3(v.x || 0, v.y || 0, v.z || 0);
        return new Vec3();
    }
    // --- Builder bindings ---

    buildSubstrate(args) {
        const { a = 2.82, nx = 10, ny = 10, nz = 3, replication = true, system, preset, step_edge } = args || {};
        
        const prev = this.activeSystemName;
        if (system) this.useSystem(system);

        let mol;
        if (preset === 'NaCl' || preset === 'NaCl(step)') {
            mol = CrystalUtils.genNaClStep({ a, nx, ny, nz });
        } else {
            // Default to NaCl step if unknown
            mol = CrystalUtils.genNaClStep({ a, nx, ny, nz });
        }

        const targetSys = this.app.systems[this.activeSystemName];
        targetSys.clear();
        for (const atom of mol.atoms) {
            targetSys.addAtom(atom.pos.x, atom.pos.y, atom.pos.z, atom.Z);
        }

        if (mol.lvec) {
            targetSys.lvec = [mol.lvec[0].clone(), mol.lvec[1].clone(), mol.lvec[2].clone()];
            const lat = this.app.getLattice(this.activeSystemName);
            lat.lvec = [mol.lvec[0].clone(), mol.lvec[1].clone(), mol.lvec[2].clone()];
            lat.nrep = { x: 1, y: 1, z: 0 };
            lat.show = true;
            this.app.updateReplicas(this.activeSystemName);
        }

        const rend = this.app.renderers[this.activeSystemName];
        if (rend) rend.update();
        this.app.requestRender();

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

    pickRandomSelection(args = {}) {
        const prev = this.activeSystemName;
        if (args.system) this.useSystem(args.system);
        const mol = this.system;
        if (!mol || !mol.selection) throw new Error('pick_random_selection: molecule or selection missing');
        const ids = Array.from(mol.selection);
        if (ids.length === 0) {
            if (args.system) this.useSystem(prev);
            throw new Error('pick_random_selection: selection is empty');
        }
        const idx = (args.pickIndex !== undefined) ? (args.pickIndex | 0) : Math.floor(Math.random() * ids.length);
        if (idx < 0 || idx >= ids.length) {
            if (args.system) this.useSystem(prev);
            throw new Error(`pick_random_selection: pickIndex out of range (0..${ids.length - 1})`);
        }
        const picked = ids[idx];
        if (args.replaceSelection) {
            mol.selection.clear();
            mol.selection.add(picked);
        }
        this.lastPickedSelection = picked;
        if (args.system) this.useSystem(prev);
        return picked;
    }

    selectBridgeCandidates(args = {}) {
        const prev = this.activeSystemName;
        if (args.system) this.useSystem(args.system);
        const mol = this.system;
        if (!mol || !mol.selectBridgeCandidates) throw new Error('select_bridge_candidates: helper not installed');
        const n = mol.selectBridgeCandidates({ requireH2: args.requireH2 });
        if (args.system) this.useSystem(prev);
        return n;
    }

    collapseBridgeAt(args = {}) {
        const prev = this.activeSystemName;
        if (args.system) this.useSystem(args.system);
        const mol = this.system;
        if (!mol || !mol.atoms) throw new Error('collapse_bridge_at: molecule missing');
        const idBridge = (args.id !== undefined) ? args.id : (args.atom !== undefined ? args.atom : this.lastPickedSelection);
        if (idBridge === undefined) {
            if (args.system) this.useSystem(prev);
            throw new Error('collapse_bridge_at: id (atom) required');
        }
        const neighIds = collapseBridgeAt(mol, idBridge);
        mol.clearSelection();
        neighIds.forEach(id => mol.selection.add(id));
        mol.dirtyExport = true;
        if (args.system) this.useSystem(prev);
    }

    collapseBridge(args = {}) {
        const prev = this.activeSystemName;
        if (args.system) this.useSystem(args.system);
        const mol = this.system;
        if (!mol || !mol.selection) throw new Error('collapse_bridge: molecule or selection missing');
        const neighIds = collapseBridgeRandom(mol);
        mol.clearSelection();
        if (neighIds) neighIds.forEach(id => mol.selection.add(id));
        mol.dirtyExport = true;
        if (args.system) this.useSystem(prev);
    }

    collapseAllBridges(args = {}) {
        const prev = this.activeSystemName;
        if (args.system) this.useSystem(args.system);
        const mol = this.system;
        if (!mol || !mol.atoms) throw new Error('collapse_all_bridges: molecule missing');
        const total = collapseAllBridges(mol, { system: args.system || this.activeSystemName });
        mol.dirtyExport = true;
        if (args.system) this.useSystem(prev);
        return total;
    }

    insertBridge(args = {}) {
        const prev = this.activeSystemName;
        if (args.system) this.useSystem(args.system);
        const mol = this.system;
        if (!mol || !mol.atoms) throw new Error('insert_bridge: molecule missing');
        const dbg = (msg) => { if (window.logger && window.logger.info) window.logger.info(msg); else console.log(msg); };
        const aId = args.aId !== undefined ? args.aId : args.a;
        const bId = args.bId !== undefined ? args.bId : args.b;
        if (aId === undefined || bId === undefined) {
            if (args.system) this.useSystem(prev);
            throw new Error('insert_bridge: aId and bId required');
        }
        const cId = insertBridge(mol, aId, bId, { ...args, dbg });
        mol.clearSelection();
        mol.selection.add(cId);
        mol.dirtyExport = true;
        if (args.system) this.useSystem(prev);
        return cId;
    }

    insertBridgeRandom(args = {}) {
        const prev = this.activeSystemName;
        if (args.system) this.useSystem(args.system);
        const mol = this.system;
        if (!mol || !mol.atoms) throw new Error('insert_bridge_random: molecule missing');
        const dbg = (msg) => { if (window.logger && window.logger.info) window.logger.info(msg); else console.log(msg); };
        const idC = insertBridgeRandom(mol, { ...args, dbg });
        mol.clearSelection();
        mol.selection.add(idC);
        mol.dirtyExport = true;
        if (args.system) this.useSystem(prev);
        return idC;
    }

}
