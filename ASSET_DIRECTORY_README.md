# Final Asset Directory

Each asset directory contains its original R, Python, and shell source chain. There is deliberately no uniform `run.py` because the assets have different upstream object contracts and terminal renderers.

Inspect a specific chain through the repository catalog:

```bash
python run_release.py --asset <asset-id>
```

The authoritative final numbering, terminal source files, output evidence, and execution stages are in `final_asset_map.tsv`.

The four final-numbered PDF supplementary-data assets retain both their analysis renderer and their recorded post-render finalization source. Those finalization steps rebuild visible content and documented transformations, not missing historical PDF bytes.
