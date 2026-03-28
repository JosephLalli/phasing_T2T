# SGDP Variation Inputs

This directory is intentionally empty in git.

The publication repository no longer tracks machine-specific symlinks to SGDP
ground-truth data. Provide these inputs either by:

1. Mounting them into the Docker container with `scripts/run_docker_smoke_test.sh`, or
2. Placing the corresponding datasets at the paths documented in `resources/README.md`.

Expected runtime paths:

- `resources/SGDP_variation/t2t`
- `resources/SGDP_variation/grch38`
