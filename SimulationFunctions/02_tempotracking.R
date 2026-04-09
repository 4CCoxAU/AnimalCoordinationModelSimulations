# ============================================================================
# Model 2: Coupled Oscillator / Tempo-Tracking
# ============================================================================
#
# Based on: Kuramoto (1984), Takahashi et al. (2013), Ermentrout (1991)
#
# Core idea: each individual has a phase oscillator, but unlike Model 1,
# the oscillators are CONTINUOUSLY COUPLED. Each individual's phase velocity
# is nudged by the partner's phase (within-cycle), and their intrinsic
# frequency gradually adapts toward the partner's frequency (across-cycle).
# No discrete resets. No refractory periods. Just smooth, mutual adjustment.
#
# Equations:
#
#   Phase dynamics (within-cycle coupling):
#     dφᵢ/dt = ωᵢ + K · sin(φⱼ − φᵢ)
#
#     K > 0:  attractive coupling → in-phase synchrony  (φⱼ − φᵢ → 0)
#     K < 0:  repulsive coupling  → anti-phase alternation (φⱼ − φᵢ → π)
#
#   Frequency adaptation (across-cycle tempo matching):
#     dωᵢ/dt = ε · (ωⱼ − ωᵢ)
#
#     ε > 0:  the slower caller speeds up, the faster caller slows down,
#             until they converge to a shared tempo.
#
#   Call emission:
#     When φᵢ crosses 2π → call emitted, φᵢ wraps to φᵢ mod 2π.
#
#   Parameters:
#     ωᵢ         angular frequency of individual i (= 2π · freq_i)
#     K          phase coupling strength
#     ε          tempo adaptation rate
#
# Key contrast with Model 1:
#   - Coupling is CONTINUOUS (every timestep), not discrete (only at calls)
#   - Both phase AND frequency are modified (not just current cycle timing)
#   - Perturbations have MEMORY: frequency changes persist across cycles
#
# Diagnostic signatures:
#   - UNIMODAL ICI (gradual tempo convergence, no mixture of interval types)
#   - Smooth phase convergence toward target (π for anti-phase)
#   - Robust entrainment across ±15-30% tempo mismatches
# ============================================================================

# Packages needed to run the script:
pacman::p_load(ggplot2, dplyr, patchwork, cowplot)

# ============================================================================
# Simulation function
# ============================================================================

simulate_coupled_oscillator <- function(
    duration = 500,            # Total simulation time (seconds)
    freq1_init = 0.8,         # Initial frequency of individual 1 (Hz)
    freq2_init = 1.2,         # Initial frequency of individual 2 (Hz)
    K = -0.1,                 # Phase coupling strength (K < 0 = anti-phase)
    epsilon = 0.05,           # Tempo adaptation rate (how fast freqs converge)
    dt = 0.01,                # Integration timestep (seconds)
    phase_noise_sd = 0.08,    # Noise in phase dynamics per timestep
    freq_noise_sd = 0.001,    # Noise in frequency adaptation per timestep
    d_call = 0.2,             # Mean call duration (seconds)
    duration_noise_sd = 0.1,  # Noise on call duration (SD)
    ici_noise_sd = 0.05,      # Timing jitter added to call onset (SD)
    initial_phase_diff = 0.3 * pi  # Starting phase offset between the two
) {

  # ---- Set up time grid ----
  # Unlike Model 1 (event-driven), this model uses fixed-timestep Euler
  # integration because the coupling is CONTINUOUS: both oscillators
  # influence each other at every moment, not just at call events.
  times   <- seq(0, duration, by = dt)
  n_steps <- length(times)

  # ---- State arrays ----
  # We track the full trajectory of phases and frequencies so we can
  # inspect convergence dynamics later if needed.
  phases <- matrix(0, nrow = n_steps, ncol = 2)  # columns = individuals
  freqs  <- matrix(0, nrow = n_steps, ncol = 2)

  # Initial conditions
  phases[1, 1] <- 0
  phases[1, 2] <- initial_phase_diff
  freqs[1, 1]  <- freq1_init
  freqs[1, 2]  <- freq2_init

  # ---- Output storage ----
  voc_times_1     <- numeric()
  voc_times_2     <- numeric()
  voc_durations_1 <- numeric()
  voc_durations_2 <- numeric()

  # ---- Main integration loop (Euler method) ----
  for (i in 2:n_steps) {

    # ---- Frequency adaptation ----
    # Each individual's frequency moves toward the partner's frequency.
    # This is the "tempo-tracking" part: the faster caller slows down,
    # the slower caller speeds up, until they match.
    freq_diff <- freqs[i-1, 2] - freqs[i-1, 1]   # how much faster is ind 2?
    freqs[i, 1] <- freqs[i-1, 1] + epsilon * freq_diff +
                   rnorm(1, 0, freq_noise_sd)
    freqs[i, 2] <- freqs[i-1, 2] - epsilon * freq_diff +
                   rnorm(1, 0, freq_noise_sd)

    # Keep frequencies positive and within a reasonable range
    freqs[i, 1] <- max(0.3, min(2.5, freqs[i, 1]))
    freqs[i, 2] <- max(0.3, min(2.5, freqs[i, 2]))

    # ---- Phase coupling (Kuramoto) ----
    # The key equation: dφᵢ/dt = ωᵢ + K · sin(φⱼ − φᵢ)
    #
    # With K < 0 (repulsive coupling):
    #   If φⱼ − φᵢ < π: sin > 0, but K < 0 → slows φᵢ → pushes apart
    #   If φⱼ − φᵢ > π: sin < 0, but K < 0 → speeds φᵢ → pulls together
    #   Net effect: phase difference converges to π (anti-phase alternation)
    #
    # With K > 0 (attractive coupling):
    #   Phase difference converges to 0 (in-phase synchrony)

    phase_diff_12 <- phases[i-1, 2] - phases[i-1, 1]  # φ₂ − φ₁
    phase_diff_21 <- phases[i-1, 1] - phases[i-1, 2]  # φ₁ − φ₂

    # Convert frequency (Hz) to angular velocity (rad/s): ω = 2π·f
    dphase_1 <- 2 * pi * freqs[i, 1] +              # intrinsic velocity
                K * sin(phase_diff_12) +              # coupling from ind 2
                rnorm(1, 0, phase_noise_sd)           # phase noise
    dphase_2 <- 2 * pi * freqs[i, 2] +              # intrinsic velocity
                K * sin(phase_diff_21) +              # coupling from ind 1
                rnorm(1, 0, phase_noise_sd)           # phase noise

    # Euler step: advance phases
    phases[i, 1] <- phases[i-1, 1] + dphase_1 * dt
    phases[i, 2] <- phases[i-1, 2] + dphase_2 * dt

    # ---- Vocalization detection ----
    # A call is emitted when the phase crosses the 2π threshold.
    # This is the equivalent of "the oscillator fires" -- but unlike
    # Model 1, the phase just wraps around; there's no discrete reset
    # or change in dynamics.

    if (phases[i, 1] >= 2 * pi && phases[i-1, 1] < 2 * pi) {
      # Individual 1's phase just crossed 2π → emit a call
      # Add small timing jitter to the onset time (biological variability)
      voc_time <- times[i] + rnorm(1, 0, ici_noise_sd)
      voc_times_1     <- c(voc_times_1, voc_time)
      voc_durations_1 <- c(voc_durations_1,
                           max(0.01, d_call + rnorm(1, 0, duration_noise_sd)))
      # Wrap phase (NOT reset -- just modular arithmetic)
      phases[i, 1] <- phases[i, 1] %% (2 * pi)
    }

    if (phases[i, 2] >= 2 * pi && phases[i-1, 2] < 2 * pi) {
      # Individual 2's phase just crossed 2π → emit a call
      voc_time <- times[i] + rnorm(1, 0, ici_noise_sd)
      voc_times_2     <- c(voc_times_2, voc_time)
      voc_durations_2 <- c(voc_durations_2,
                           max(0.01, d_call + rnorm(1, 0, duration_noise_sd)))
      phases[i, 2] <- phases[i, 2] %% (2 * pi)
    }

    # Wrap phases to [0, 2π) regardless (handles accumulation over time)
    phases[i, 1] <- phases[i, 1] %% (2 * pi)
    phases[i, 2] <- phases[i, 2] %% (2 * pi)
  }

  # ---- Return results ----
  list(
    voc_times_1     = voc_times_1,
    voc_times_2     = voc_times_2,
    voc_durations_1 = voc_durations_1,
    voc_durations_2 = voc_durations_2,
    # Also return full trajectories for optional deeper inspection
    times  = times,
    phases = phases,   # n_steps x 2 matrix
    freqs  = freqs,    # n_steps x 2 matrix
    params = list(
      freq1_init = freq1_init, freq2_init = freq2_init,
      K = K, epsilon = epsilon, dt = dt,
      phase_noise_sd = phase_noise_sd, freq_noise_sd = freq_noise_sd,
      d_call = d_call, duration_noise_sd = duration_noise_sd,
      ici_noise_sd = ici_noise_sd,
      initial_phase_diff = initial_phase_diff,
      duration = duration
    )
  )
}
