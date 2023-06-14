import pytest


@pytest.mark.skip(
    reason="Application is a daemon process, does not return anything to test."
)
def test():
    from aware.app import Application

    app = Application()
    app.run()


if __name__ == "__main__":
    test()
