import pathlib

import pytest
import tomlkit as tomllib
from click.testing import CliRunner

from plume.cli import plume


def test_command_line_interface():
    runner = CliRunner()
    result = runner.invoke(plume, ["--help"])
    assert result.exit_code == 0

    result = runner.invoke(plume, ["--version"])
    assert result.exit_code == 0
    assert "version" in result.output


@pytest.mark.parametrize("subcommand", ("generate", "run", "setup"))
def test_subcommand_help(subcommand):
    runner = CliRunner()
    result = runner.invoke(plume, [subcommand, "--help"])
    assert result.exit_code == 0


def test_setup(tmpdir):
    runner = CliRunner()
    with tmpdir.as_cwd():
        result = runner.invoke(plume, ["setup"])
        assert result.exit_code == 0
        assert pathlib.Path("plume.toml").is_file()

        with open("plume.toml") as fp:
            params = tomllib.load(fp)
        assert "plume" in params
