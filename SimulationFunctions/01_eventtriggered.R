# ============================================================================
# Model 1: Inhibitory Reset (Event-Triggered)
# ============================================================================
#
# Based on: Greenfield (1994), Sismondo (1990)
#
# Core idea: each individual has a free-running oscillator that cycles
# from phase 0 to 2π. When it hits 2π, a call is emitted. When a partner
# calls, the oscillator is RESET back to 0 (inhibitory resetting), but
# only if the individual is past a refractory window. After reset, the
# next cycle is shortened (post-inhibitory rebound).
#
# Equations:
#
#   Free-running:  dφᵢ/dt = 2π/Tᵢ        (phase advances at constant rate)
#   Call emission:  φᵢ = 2π → call, φᵢ → 0 (call fires at threshold)
#   Partner call arrives when focal is at phase φᵢ:
#     If φᵢ > φ_ref:  φᵢ → 0, next period = d + ρTᵢ  (reset + rebound)
#     If φᵢ ≤ φ_ref:  no response (refractory)
#
#   Parameters:
#     Tᵢ           natural period of individual i
#     ρ ∈ (0,1)    rebound ratio (how much the period shortens after reset)
#     d            call duration (also sets the inhibition delay)
#     φ_ref        refractory threshold = 2π(d/Tᵢ), i.e., the first d/T
#                  fraction of the cycle is refractory to stimuli
#
# Diagnostic signatures:
#   - Bimodal ICI: mode at T (unperturbed) and ~2(d + ρT) (alternation)
#   - Oscillatory phase jitter (not smooth convergence)
#   - Perturbations affect only the current cycle; intrinsic ω never changes
# ============================================================================

# Packages needed to run the script:
pacman::p_load(ggplot2, dplyr, patchwork, cowplot)

# ============================================================================
# Simulation function
# ============================================================================

simulate_event_triggered <- function(
    duration = 500,          # Total simulation time (seconds)
    T1 = 1.0,               # Natural period of individual 1 (seconds)
    T2 = 1.0,               # Natural period of individual 2 (seconds)
    rho = 0.7,              # Rebound ratio: after reset, period = d + rho*T
    d_call = 0.15,           # Call duration / inhibition delay (seconds)
    noise_sd = 0.2,         # Gaussian noise added to each period (SD)
    duration_noise_sd = 0.1 # Gaussian noise on call duration (SD)
) {

  # ---- Refractory threshold ----
  # φ_ref = fraction of cycle occupied by the call, mapped to radians.
  # If the partner's call arrives while the focal individual is still in
  # this early part of the cycle (just called), the stimulus is ignored.
  phi_ref_1 <- 2 * pi * (d_call / T1)
  phi_ref_2 <- 2 * pi * (d_call / T2)

  # ---- Output storage ----
  # We store the onset time and duration of every call produced.
  voc_times_1     <- numeric()
  voc_times_2     <- numeric()
  voc_durations_1 <- numeric()
  voc_durations_2 <- numeric()

  # ---- State variables ----
  time    <- 0                              # Current simulation clock
  phase_1 <- runif(1, 0, 2 * pi)           # Initial phase (random start)
  phase_2 <- runif(1, 0, 2 * pi)           # Initial phase (random start)
  period_1 <- max(0.1, T1 + rnorm(1, 0, noise_sd))  # Current cycle length
  period_2 <- max(0.1, T2 + rnorm(1, 0, noise_sd))  # Current cycle length

  # ---- Main event-driven loop ----
  # At each step we ask: "who reaches phase 2π first?"
  # That individual fires a call. Then we check whether the partner
  # should be reset.

  while (time < duration) {

    # How long until each individual's phase reaches 2π?
    # phase advances at rate 2π/period, so remaining time =
    # (2π - current_phase) / (2π / period) = period * (1 - phase/2π)
    time_to_call_1 <- period_1 * (1 - phase_1 / (2 * pi))
    time_to_call_2 <- period_2 * (1 - phase_2 / (2 * pi))

    if (time_to_call_1 <= time_to_call_2) {

      # ============================================================
      # Individual 1 calls first
      # ============================================================

      # Advance the clock by the time until individual 1's call
      dt   <- time_to_call_1
      time <- time + dt
      if (time >= duration) break

      # Record the call (onset time + noisy duration)
      dur <- max(0.01, d_call + rnorm(1, 0, duration_noise_sd))
      voc_times_1     <- c(voc_times_1, time)
      voc_durations_1 <- c(voc_durations_1, dur)

      # ---- Effect on individual 2 ----
      # First, advance individual 2's phase by the elapsed time dt.
      # They've been free-running while waiting for individual 1 to call.
      phase_2 <- phase_2 + dt * (2 * pi / period_2)

      # Now check the refractory condition:
      # Is individual 2 past the refractory window?
      if (phase_2 > phi_ref_2) {
        # YES: individual 2 gets RESET.
        # Their phase goes back to 0 and their next period is shortened
        # to d_call (must wait for the inhibiting call to finish)
        # plus rho * T2 (the post-inhibitory rebound period).
        phase_2  <- 0
        period_2 <- max(0.1, d_call + rho * T2 + rnorm(1, 0, noise_sd))
      }
      # If phase_2 <= phi_ref_2: individual 2 is refractory (just called
      # recently), so the stimulus has no effect. Phase and period unchanged.

      # ---- Individual 1 starts a new normal cycle ----
      # The calling individual always returns to their natural period T.
      # This is key: perturbations affect only the PARTNER's current cycle,
      # never the caller's intrinsic frequency.
      phase_1  <- 0
      period_1 <- max(0.1, T1 + rnorm(1, 0, noise_sd))

    } else {

      # ============================================================
      # Individual 2 calls first (mirror logic)
      # ============================================================

      dt   <- time_to_call_2
      time <- time + dt
      if (time >= duration) break

      dur <- max(0.01, d_call + rnorm(1, 0, duration_noise_sd))
      voc_times_2     <- c(voc_times_2, time)
      voc_durations_2 <- c(voc_durations_2, dur)

      # Advance individual 1's phase by dt
      phase_1 <- phase_1 + dt * (2 * pi / period_1)

      # Refractory check on individual 1
      if (phase_1 > phi_ref_1) {
        # RESET individual 1
        phase_1  <- 0
        period_1 <- max(0.1, d_call + rho * T1 + rnorm(1, 0, noise_sd))
      }

      # Individual 2 starts a new normal cycle
      phase_2  <- 0
      period_2 <- max(0.1, T2 + rnorm(1, 0, noise_sd))
    }
  }

  # ---- Return results ----
  list(
    voc_times_1     = voc_times_1,
    voc_times_2     = voc_times_2,
    voc_durations_1 = voc_durations_1,
    voc_durations_2 = voc_durations_2,
    params = list(
      T1 = T1, T2 = T2, rho = rho, d_call = d_call,
      phi_ref_1 = phi_ref_1, phi_ref_2 = phi_ref_2,
      noise_sd = noise_sd, duration_noise_sd = duration_noise_sd,
      duration = duration
    )
  )
}
