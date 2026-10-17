"""
wormcat3.cli.main
~~~~~~~~~~~~~~~~~
Top-level CLI execution entrypoint.
"""

import logging
import sys
from pathlib import Path
from typing import Optional, Sequence

from pydantic import ValidationError

from wormcat3.cli.commands import (
    handle_annotations_list,
    handle_batch,
    handle_enrich,
    handle_extract,
    handle_gmt,
    handle_gsea,
)
from wormcat3.cli.models import (
    BatchConfig,
    EnrichmentConfig,
    ExtractSheetsConfig,
    GMTExportConfig,
    GSEAConfig,
)
from wormcat3.cli.parser import build_parser
from wormcat3.logger import configure_logging, get_logger
from wormcat3.wormcat_error import WormcatError

logger = get_logger(__name__)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entrypoint invoked from console script or direct call."""
    parser = build_parser()
    parsed_args = parser.parse_args(args=argv)

    # Configure logging level based on global flags
    log_level = logging.INFO
    if parsed_args.verbose:
        log_level = logging.DEBUG
    elif parsed_args.quiet:
        log_level = logging.WARNING

    configure_logging(level=log_level, log_file=parsed_args.log_file)

    try:
        if parsed_args.subcommand == "enrich":
            enrich_cfg = EnrichmentConfig(
                input_file=Path(parsed_args.input),
                background_file=Path(parsed_args.background) if parsed_args.background else None,
                annotation_file=parsed_args.annotation,
                method=parsed_args.method,
                threshold=parsed_args.threshold,
                gene_type=parsed_args.gene_type,
                output_dir=Path(parsed_args.output_dir),
                title=parsed_args.title,
                email=parsed_args.email,
                no_viz=parsed_args.no_viz,
                verbose=parsed_args.verbose,
                quiet=parsed_args.quiet,
                log_file=Path(parsed_args.log_file) if parsed_args.log_file else None,
            )
            return handle_enrich(enrich_cfg)

        elif parsed_args.subcommand == "batch":
            batch_cfg = BatchConfig(
                input_path=Path(parsed_args.input),
                background_file=Path(parsed_args.background) if parsed_args.background else None,
                annotation_file=parsed_args.annotation,
                method=parsed_args.method,
                threshold=parsed_args.threshold,
                gene_type=parsed_args.gene_type,
                output_dir=Path(parsed_args.output_dir),
                title=parsed_args.title,
                email=parsed_args.email,
                verbose=parsed_args.verbose,
                quiet=parsed_args.quiet,
                log_file=Path(parsed_args.log_file) if parsed_args.log_file else None,
            )
            return handle_batch(batch_cfg)

        elif parsed_args.subcommand == "gsea":
            gsea_cfg = GSEAConfig(
                input_file=Path(parsed_args.input),
                annotation_file=parsed_args.annotation,
                output_dir=Path(parsed_args.output_dir),
                title=parsed_args.title,
                email=parsed_args.email,
                min_size=parsed_args.min_size,
                max_size=parsed_args.max_size,
                permutation_num=parsed_args.permutations,
                weight=parsed_args.weight,
                seed=parsed_args.seed,
                threads=parsed_args.threads,
                verbose=parsed_args.verbose,
                quiet=parsed_args.quiet,
                log_file=Path(parsed_args.log_file) if parsed_args.log_file else None,
            )
            return handle_gsea(gsea_cfg)

        elif parsed_args.subcommand == "gmt":
            gmt_cfg = GMTExportConfig(
                output_dir=Path(parsed_args.output_dir),
                annotation_file=parsed_args.annotation,
                gene_type=parsed_args.gene_type,
                prefix=parsed_args.prefix,
            )
            return handle_gmt(gmt_cfg)

        elif parsed_args.subcommand == "annotations":
            if parsed_args.annot_action == "list" or parsed_args.annot_action is None:
                return handle_annotations_list()
            else:
                parser.error(f"Unknown annotations action: {parsed_args.annot_action}")

        elif parsed_args.subcommand == "extract":
            extract_cfg = ExtractSheetsConfig(
                input_file=Path(parsed_args.input),
                output_dir=Path(parsed_args.output_dir),
            )
            return handle_extract(extract_cfg)

        else:
            parser.print_help(sys.stderr)
            return 1

    except ValidationError as err:
        logger.error(f"Configuration Validation Error:\n{err}")
        return 1
    except WormcatError as err:
        logger.error(f"WormCat Error: {err}")
        return 1
    except Exception as err:
        if parsed_args.verbose:
            logger.exception(f"Unexpected fatal error: {err}")
        else:
            logger.error(f"Fatal error: {err} (use -v/--verbose for full traceback)")
        return 1


if __name__ == "__main__":
    sys.exit(main())
