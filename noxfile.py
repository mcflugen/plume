import pathlib
import shutil
from itertools import chain

import nox

ROOT = pathlib.Path(__file__).parent


@nox.session(venv_backend="mamba")
def install(session: nox.Session) -> None:
    """Install the package."""
    session.conda_install(
        "--file=requirements.txt",
        "gsl",
        channel=["nodefaults", "conda-forge"],
    )
    session.install("-e", ".", "--no-deps")


@nox.session(venv_backend="mamba")
def test(session: nox.Session) -> None:
    """Run the tests."""
    install(session)
    session.conda_install(
        "--file=requirements-testing.txt", channel=["nodefaults", "conda-forge"]
    )

    session.run("pytest", "--cov=src/plume", "-vvv")
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

    install(session)
    session.conda_install(
        "notebook",
        "--file=requirements-testing.txt",
        channel=["nodefaults", "conda-forge"],
    )
    session.install("git+https://github.com/mcflugen/nbmake.git@mcflugen/add-markers")

    session.run(*args)


@nox.session(name="test-cli", venv_backend="mamba")
def test_cli(session: nox.Session) -> None:
    """Test the command line interface."""
    install(session)

    session.run("plume", "--version")
    session.run("plume", "--help")
    session.run("plume", "generate", "--help")
    session.run("plume", "setup", "--help")
    session.run("plume", "run", "--help")
    session.run("plume", "generate", "plume.toml")


@nox.session
def lint(session: nox.Session) -> None:
    """Look for lint."""
    skip_hooks = [] if "--no-skip" in session.posargs else ["check-manifest", "pyroma"]

    session.install("pre-commit")
    session.run("pre-commit", "run", "--all-files", env={"SKIP": ",".join(skip_hooks)})


@nox.session
def towncrier(session: nox.Session) -> None:
    """Check that there is a news fragment."""
    session.install("towncrier")
    session.run("towncrier", "check", "--compare-with", "origin/main")


@nox.session(name="sync-requirements", python="3.11", venv_backend="conda")
def sync_requirements(session: nox.Session) -> None:
    """Sync requirements.in with pyproject.toml."""
    with open("requirements.txt", "w") as fp:
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
    PROJECT = "sed_plume"
    ROOT = pathlib.Path(__file__).parent

    shutil.rmtree("build", ignore_errors=True)
    shutil.rmtree("wheelhouse", ignore_errors=True)
    shutil.rmtree(f"src/{PROJECT}.egg-info", ignore_errors=True)
    shutil.rmtree(".pytest_cache", ignore_errors=True)
    shutil.rmtree(".venv", ignore_errors=True)
    for p in chain(
        ROOT.rglob("*.c"),
        ROOT.rglob("*.py[co]"),
        ROOT.rglob("*.so"),
        ROOT.rglob("__pycache__"),
    ):
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
