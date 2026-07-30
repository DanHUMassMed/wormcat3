import wormcat3


def test_version():
    assert hasattr(wormcat3, "__version__")
    assert isinstance(wormcat3.__version__, str)


def test_imports():
    assert hasattr(wormcat3, "Wormcat")
    assert hasattr(wormcat3, "AnnotationsManager")
    assert hasattr(wormcat3, "PAdjustMethod")
    assert hasattr(wormcat3, "WormcatError")
