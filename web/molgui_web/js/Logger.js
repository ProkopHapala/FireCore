class Logger {
    constructor() {
        this.domElement = null; // Will be set by GUI
        this.level = 1; // 0: DEBUG, 1: INFO, 2: WARN, 3: ERROR
        this.levels = { 'DEBUG': 0, 'INFO': 1, 'WARN': 2, 'ERROR': 3 };
    }

    setContainer(element) {
        this.domElement = element;
    }

    setLevel(levelStr) {
        if (this.levels.hasOwnProperty(levelStr)) {
            this.level = this.levels[levelStr];
            this.info(`Log Level set to ${levelStr}`);
        }
    }

    clear() {
        if (this.domElement) {
            this.domElement.textContent = '';
        }
    }

    log(message, level = 'INFO') {
        const lvl = this.levels[level] !== undefined ? this.levels[level] : 1;
        if (lvl < this.level) return;

        const timestamp = new Date().toLocaleTimeString();
        const formattedMsg = `[${timestamp}] [${level}] ${message}`;

        // Console
        if (level === 'ERROR') console.error(formattedMsg);
        else if (level === 'WARN') console.warn(formattedMsg);
        else console.log(formattedMsg);

        // DOM
        if (this.domElement) {
            const line = document.createElement('div');
            line.textContent = formattedMsg;
            if (level === 'ERROR') line.style.color = '#ff5555';
            else if (level === 'WARN') line.style.color = '#ffaa00';
            else if (level === 'DEBUG') line.style.color = '#55ff55';

            this.domElement.appendChild(line);
            this.domElement.scrollTop = this.domElement.scrollHeight;
        }
    }

    debug(msg) { this.log(msg, 'DEBUG'); }
    info(msg) { this.log(msg, 'INFO'); }
    warn(msg) { this.log(msg, 'WARN'); }
    error(msg) { this.log(msg, 'ERROR'); }
}

// Global instance
window.logger = new Logger();
