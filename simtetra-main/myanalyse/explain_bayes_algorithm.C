Initial guess: π(j)

LOOP:
    For each measured bin i:
        Compute P(j|i) = R_ij * π(j) / sum_k R_ik * π(k)

    For each true bin j:
        f_new(j) = Σ_i g(i) * P(j|i)

    Normalize, enforce positivity

    π = f_new   (becomes prior)

    if (χ² changes slowly):
          STOP