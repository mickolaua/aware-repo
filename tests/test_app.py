from aware.app import Application
import pytest


@pytest.mark.skip("Application is a daemon process, does not return anything to test.")
def test():
    app = Application()
    app.run()


if __name__ == "__main__":
    test()
