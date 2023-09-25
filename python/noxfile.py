import nox


@nox.session
def tests(session):
    session.install("-e", "..[test]")
    session.run("pytest", "-v", "test")
