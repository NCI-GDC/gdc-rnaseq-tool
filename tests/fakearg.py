from typing import Any


class FakeArgs:
    """Mock args"""

    def __setattr__(self, key: str, value: Any) -> None:
        self.__dict__[key] = value
