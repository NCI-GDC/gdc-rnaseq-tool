
class FakeArgs:
    """Mock args"""
    def __setattr__(self, key, value):
        self.__dict__[key] = value
