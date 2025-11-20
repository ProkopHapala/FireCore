class Logger {
    static LEVEL_DEBUG = 0;
    static LEVEL_INFO = 1;
    static LEVEL_WARN = 2;
    static LEVEL_ERROR = 3;

    constructor() {
        this.domElement = document.getElementById('log-content');
        this.level = Logger.LEVEL_DEBUG;
    }

    log(level, message) {
        if (level < this.level) return;

        const timestamp = new Date().toLocaleTimeString();
        const msg = `[${timestamp}] ${message}`;

        // Console output
        switch(level) {
            case Logger.LEVEL_ERROR: console.error(msg); break;
            case Logger.LEVEL_WARN: console.warn(msg); break;
            case Logger.LEVEL_INFO: console.info(msg); break;
            default: console.log(msg);
        }

        // DOM output
        if (this.domElement) {
            const line = document.createElement('div');
            line.textContent = msg;
            
            switch(level) {
                case Logger.LEVEL_ERROR: line.className = 'log-error'; break;
                case Logger.LEVEL_WARN: line.className = 'log-warn'; break;
                case Logger.LEVEL_INFO: line.className = 'log-info'; break;
                default: line.className = 'log-debug';
            }
            
            this.domElement.appendChild(line);
            this.domElement.scrollTop = this.domElement.scrollHeight;
        }
    }

    debug(msg) { this.log(Logger.LEVEL_DEBUG, msg); }
    info(msg) { this.log(Logger.LEVEL_INFO, msg); }
    warn(msg) { this.log(Logger.LEVEL_WARN, msg); }
    error(msg) { this.log(Logger.LEVEL_ERROR, msg); }
}

// Global instance
window.logger = new Logger();
