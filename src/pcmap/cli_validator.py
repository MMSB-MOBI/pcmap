import click
import pypstruct
from pypstruct.coordinates import Structure
from rich import print

from .io import tripletParser


def parse_PDB_file(ctx, opts, value: str) -> Structure:
    _ = pypstruct.parseFilePDB(value)
    if _ is None:
        raise click.BadParameter(f"Could not parse PDB file '{value}'")
    return _


def positive_float(ctx, param, value):
    try:
        fvalue = float(value)
    except ValueError:
        raise click.BadParameter("Value must be a valid float.")
    if fvalue <= 0:
        raise click.BadParameter("Value must be a positive float.")
    return fvalue


def positive_int(ctx, param, value):
    try:
        fvalue = int(value)
    except ValueError:
        raise click.BadParameter("Value must be a valid integer.")
    if fvalue <= 0:
        raise click.BadParameter("Value must be a positive integer.")
    return fvalue


def check_rich_and_atomic(ctx, param, value):
    if ctx.params.get("rich"):
        if not ctx.params.get("atomic"):
            raise click.UsageError(
                "Rich contact map output (--rich) requires the atomic flag (--atomic)"
            )
    return value


def parse_triplet(ctx, param, value) -> None | tuple[float, float, float]:
    if value is None:
        return value
    try:
        _ = tripletParser(value)
        return _
    except:
        raise click.BadParameter(f"{param} value must be float triplet")


import os


def path_writable(ctx, param, file_path):
    # If the file does not exist, check if the directory is writable
    if not os.path.exists(file_path):
        dirname = os.path.dirname(file_path)
        # If directory does not exist, return False
        if not os.path.exists(dirname):
            raise click.BadParameter(f"Output folder doesnt exist '{dirname}'")
        if not os.access(dirname, os.W_OK):
            raise click.BadParameter(f"Output folder is not writable '{dirname}'")
    # If the file exists, check if it is writable
    if not os.access(file_path, os.W_OK):
        raise click.BadParameter(f"Output file is invalid '{file_path}'")

    return file_path
