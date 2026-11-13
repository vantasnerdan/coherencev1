# coherencev1

A coarse end-to-end model of the coherence-domain hypothesis for hydrogen-loaded metals (Ni/H, Pd/D), refined in passes.

The aim is a model that maps apparatus-level inputs to measured outputs, so there is something to engineer against. It does not set out to show the effect exists.

## Method

Model the whole chain in its cheapest form first, then add resolution where the end-to-end output is most sensitive. Passes are sized by compute.

- **Pass 0, negligible compute:** a minimal domain model (matter coupled to one field mode, a few ODEs), looked-up material scales, and one-line placeholders for the rest of the chain.
- **Pass 1, light compute:** sweep the pass 0 model, compress it into closed-form laws (PySR, SINDy), rank the unknowns by sensitivity.
- **Pass 2, moderate compute:** one number from a real electronic-structure run for the top-ranked unknown; refit.
- **Pass 3 and on, heavy compute:** converged runs where the ranking points.

## Layout

One directory per work item, named after it (`p0_scales/`, `p0_domain/`, ...). Each holds a small self-contained script and its output table.
