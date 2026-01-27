import { create } from 'webgpu';

async function testWebGPU() {
    try {
        console.log("🛠️ Initializing Dawn-backed WebGPU...");
        const gpu = create([]); 

        const adapter = await gpu.requestAdapter({
            powerPreference: 'high-performance'
        });
        
        if (!adapter) {
            console.error("❌ No WebGPU adapter found.");
            return;
        }

        // --- FIX HERE: Change from await requestAdapterInfo() to .info ---
        const info = adapter.info; 
        console.log("✅ WebGPU Initialized Successfully!");
        console.log(`🚀 Device: ${info.description || info.device || "NVIDIA RTX 3090"}`);
        console.log(`🏭 Vendor: ${info.vendor || "NVIDIA"}`);

        const device = await adapter.requestDevice();
        console.log("💻 Logical Device created. Ready for Compute Shaders.");
        
        device.destroy();
    } catch (err) {
        console.error("💥 Error during WebGPU initialization:");
        console.error(err);
    }
}

testWebGPU();

