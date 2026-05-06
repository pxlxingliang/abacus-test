import sys

_LAZY_IMPORTS = {
    'ReadInput': 'abacustest.lib_prepare.abacus',
    'WriteInput': 'abacustest.lib_prepare.abacus',
    'AbacusStru': 'abacustest.lib_prepare.abacus',
    'ReadKpt': 'abacustest.lib_prepare.abacus',
    'WriteKpt': 'abacustest.lib_prepare.abacus',
    'RESULT': 'abacustest.lib_collectdata.collectdata',
    'AbacusSTRU': 'abacustest.lib_prepare.stru',
}

def __getattr__(name):
    if name in _LAZY_IMPORTS:
        import importlib
        module = importlib.import_module(_LAZY_IMPORTS[name])
        globals()[name] = getattr(module, name)
        return globals()[name]
    raise AttributeError(name)

__all__ = list(_LAZY_IMPORTS.keys())