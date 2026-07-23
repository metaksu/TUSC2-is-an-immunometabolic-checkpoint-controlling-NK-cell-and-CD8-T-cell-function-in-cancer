"""Import guard: every module in the package imports cleanly, and every symbol a
module imports from its own config actually exists in that config.

The config check uses a static AST parse, so it also catches names imported
inside a function body (which a plain import-the-module smoke check would miss),
while excluding path imports that legitimately come from elsewhere.
"""
import ast
import importlib
from pathlib import Path

import pytest

PKG = "tusc2_deg.nk"
PKG_DIR = Path(importlib.import_module(PKG).__file__).resolve().parent


def _package_modules() -> list[str]:
    mods = []
    for p in sorted(PKG_DIR.glob("*.py")):
        if p.name == "__init__.py":
            continue
        mods.append(f"{PKG}.{p.stem}")
    return mods


@pytest.mark.parametrize("modname", _package_modules())
def test_module_imports_cleanly(modname):
    """Every .py in the package imports without error (catches module-level
    breakage, including the heavy sensitivity script's top-level imports)."""
    importlib.import_module(modname)


_PKG_CONFIG_MODULES = {f"{PKG}.config", ".config", "config"}


def _pkg_config_imports_in(path: Path) -> set[str]:
    """Names imported specifically from the PACKAGE config
    (`from tusc2_deg.nk.config import ...` or `from .config`),
    anywhere in the file (module-level OR function-local), via static AST parse —
    so we never execute the heavy main() body.

    Imports from the data-paths module (`from tusc2_deg.paths import ...`) are
    EXCLUDED: those names (paths) legitimately do not live in the package config.
    We disambiguate by node.level (relative `.config` => level 1) and by the
    fully-qualified package path.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    names: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.ImportFrom):
            continue
        mod = node.module or ""
        is_pkg_config = (
            mod == f"{PKG}.config"                      # absolute package config
            or (node.level >= 1 and mod == "config")    # relative .config
        )
        if is_pkg_config:
            for alias in node.names:
                if alias.name != "*":
                    names.add(alias.name)
    return names


@pytest.mark.parametrize("modname", _package_modules())
def test_config_symbols_resolve(modname):
    """Every symbol a package module imports from ITS OWN (package) config must
    EXIST in that config — even when the import is buried inside a function body.
    Root-config path imports are excluded by _pkg_config_imports_in."""
    pkg_config = importlib.import_module(f"{PKG}.config")
    path = PKG_DIR / (modname.split(".")[-1] + ".py")
    for name in _pkg_config_imports_in(path):
        assert hasattr(pkg_config, name), \
            f"{modname} imports {PKG}.config.{name} which does not exist"
