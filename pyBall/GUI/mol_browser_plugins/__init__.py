"""Built-in VispyMolBrowser plugins."""
from pyBall.GUI.mol_browser_plugins.base import MolBrowserContext, MolBrowserPlugin
from pyBall.GUI.mol_browser_plugins.registry import MolBrowserPluginHost, MolBrowserPluginRegistry
from pyBall.GUI.mol_browser_plugins.vibration_spectrum import VibrationSpectrumPlugin


def default_plugin_registry() -> MolBrowserPluginRegistry:
    reg = MolBrowserPluginRegistry()
    reg.register(VibrationSpectrumPlugin())
    return reg

__all__ = [
    'MolBrowserContext',
    'MolBrowserPlugin',
    'MolBrowserPluginHost',
    'MolBrowserPluginRegistry',
    'VibrationSpectrumPlugin',
    'default_plugin_registry',
]
