from __future__ import annotations

# from aware import plugin_manager
from aware.alert import AlertParsers


def test():
    plugins = AlertParsers.values()
    assert len(plugins) > 1, "Alert parsers not loaded!"
    
    
if __name__ == "__main__":
    test()

