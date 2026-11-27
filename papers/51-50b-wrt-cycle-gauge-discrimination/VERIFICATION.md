# Paper 51+50B verification

- The cycle asymptotic covers each fixed even cycle and both shifted-level parities. The separate finite-level modular-matrix reproduction is in the [derivation supplement](supplement/README.txt).
- The 24-generator cap–fiber presentation and the eight-generator reduction are independently checked by [exact Smith-form replay](check_signed_homology.py).
- The order-21 result uses exact field arithmetic. [Replay](check_order21.py) compares two rank reductions for each character; see [frozen counts](order21-replay.txt). The proof and scope were cross-reviewed internally.
- The [paper PDF](manuscript.pdf) and [source](main.tex) are recorded in [SHA256SUMS](SHA256SUMS). There has been no external peer review.
