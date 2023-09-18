import pathlib
import shutil
from itertools import chain

import nox

ROOT = pathlib.Path(__file__).parent


@nox.session(venv_backend="mamba")
def test(session: nox.Session) -> None:
    """Run the tests."""
    session.conda_install(
        "--file=requirements.txt", "--file=requirements-testing.txt", "gsl"
    )
    session.install(".", "--no-deps")

    session.run("pytest", "--cov=plume", "-vvv")
    session.run("coverage", "report", "--ignore-errors", "--show-missing")


@nox.session(name="test-notebooks", venv_backend="mamba")
def test_notebooks(session: nox.Session) -> None:
    """Run the notebooks."""
    args = [
        "pytest",
        "notebooks",
        "--nbmake",
        "--nbmake-kernel=python3",
        "--nbmake-timeout=3000",
        "-n",
        "auto",
        "-vvv",
    ] + session.posargs

    session.conda_install(
        "gsl",
        "notebook",
        "--file=requirements-testing.txt",
        "--file=requirements.txt",
        channel=["nodefaults", "conda-forge"],
    )
    session.install("git+https://github.com/mcflugen/nbmake.git@mcflugen/add-markers")
    session.install("-e", ".", "--no-deps")

    session.run(*args)


@nox.session
def cli(session: nox.Session) -> None:
    """Test the command line interface."""
    session.install(".")

    session.run("plume", "--version")
    session.run("plume", "--help")
    session.run("plume", "generate", "--help")
    session.run("plume", "setup", "--help")
    session.run("plume", "run", "--help")
    session.run("plume", "generate", "plume.toml")


@nox.session
def lint(session: nox.Session) -> None:
    """Look for lint."""
    session.install("pre-commit")
    session.run("pre-commit", "run", "--all-files")


@nox.session
def towncrier(session: nox.Session) -> None:
    """Check that there is a news fragment."""
    session.install("towncrier")
    session.run("towncrier", "check", "--compare-with", "origin/main")


@nox.session(name="sync-requirements", python="3.11", venv_backend="conda")
def sync_requirements(session: nox.Session) -> None:
    """Sync requirements.in with pyproject.toml."""
    pypi_mapping = {"ecmwflibs": "findlibs"}

    with open("requirements.in", "w") as fp:
        session.run(
            "python",
            "-c",
            """
import os, tomllib
with open("pyproject.toml", "rb") as fp:
    print(os.linesep.join(sorted(tomllib.load(fp)["project"]["dependencies"])))
""",
            stdout=fp,
        )

    with open("requirements.in") as fp:
        pypi_requirements = set(fp.read().splitlines())

    with open("requirements-conda.in", "w") as fp:
        for requirement in sorted(
            {pypi_mapping.get(req, req) for req in pypi_requirements}
        ):
            print(requirement, file=fp)


@nox.session
def build(session: nox.Session) -> None:
    """Build sdist and wheel dists."""
    session.install("pip")
    session.install("wheel")
    session.install("setuptools")
    session.run("python", "--version")
    session.run("pip", "--version")
    session.run(
        "python", "setup.py", "bdist_wheel", "sdist", "--dist-dir", "./wheelhouse"
    )


@nox.session
def release(session):
    """Tag, build and publish a new release to PyPI."""
    session.install("zest.releaser[recommended]")
    session.install("zestreleaser.towncrier")
    session.run("fullrelease")


@nox.session
def publish_testpypi(session):
    """Publish wheelhouse/* to TestPyPI."""
    session.run("twine", "check", "wheelhouse/*")
    session.run(
        "twine",
        "upload",
        "--skip-existing",
        "--repository-url",
        "https://test.pypi.org/legacy/",
        "wheelhouse/*.tar.gz",
    )


@nox.session
def publish_pypi(session):
    """Publish wheelhouse/* to PyPI."""
    session.run("twine", "check", "wheelhouse/*")
    session.run(
        "twine",
        "upload",
        "--skip-existing",
        "wheelhouse/*.tar.gz",
    )


@nox.session(python=False)
def clean(session):
    """Remove all .venv's, build files and caches in the directory."""
    PROJECT = "plume"
    ROOT = pathlib.Path(__file__).parent

    shutil.rmtree("build", ignore_errors=True)
    shutil.rmtree("wheelhouse", ignore_errors=True)
    shutil.rmtree(f"{PROJECT}.egg-info", ignore_errors=True)
    shutil.rmtree(".pytest_cache", ignore_errors=True)
    shutil.rmtree(".venv", ignore_errors=True)
    for p in chain(ROOT.rglob("*.py[co]"), ROOT.rglob("__pycache__")):
        if p.is_dir():
            p.rmdir()
        else:
            p.unlink()


@nox.session(python=False, name="clean-docs")
def clean_docs(session: nox.Session) -> None:
    """Clean up the docs folder."""
    if (ROOT / "build" / "html").exists():
        with session.chdir(ROOT / "build"):
            shutil.rmtree("html")
