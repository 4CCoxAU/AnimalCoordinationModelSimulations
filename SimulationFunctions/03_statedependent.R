# ============================================================================
# Model 3: State-Dependent / Hidden Markov Model
# ============================================================================
#
# Based on: Popov et al. (2017), Donnarumma et al. (2017), Stowell et al. (2016)
#
# Core idea: coordination is NOT a timing problem but a DECISION problem.
# Each individual cycles through unobservable behavioral states (Active,
# Listening, Refractory) that determine calling probability. The partner's
# vocal activity modulates STATE TRANSITIONS, not phase or frequency.
# No oscillators. No intrinsic periodicity. Coordination emerges from
# probabilistic gating rather than phase-locking.
#
# States:
#   Active     – high calling probability (the animal is "ready to call")
#   Listening  – low calling probability (suppressed by partner activity)
#   Refractory – near-zero calling probability (recovery after own call)
#
# Dynamics:
#   At each timestep dt:
#     1. Transition between states based on:
#        - Current state
#        - Whether the partner called recently (within the last few timesteps)
#        - Time spent in current state (for Refractory exit)
#
#     2. Emit a call with state-dependent probability:
#        P(call | Active)     = λ_active · dt       (high)
#        P(call | Listening)  = λ_listening · dt     (low)
#        P(call | Refractory) = λ_refractory · dt    (near-zero)
#
#     3. After calling → enter Refractory (mandatory cooldown)
#        After partner calls → high probability of Active → Listening
#
# Key equations:
#
#   State transition (no partner call):
#     P(s_t | s_{t-1}) = T_baseline[s_{t-1}, s_t]
#
#   State transition (partner called):
#     P(s_t | s_{t-1}) = T_partner[s_{t-1}, s_t]
#
#   Call emission:
#     P(call | s_t) = λ_{s_t} · dt
#
#   After own call:
#     s_t → Refractory (deterministic)
#
# Parameters:
#   λ_active, λ_listening, λ_refractory  – calling rates (Hz) per state
#   T_baseline  – 3×3 transition matrix (no partner call)
#   T_partner   – 3×3 transition matrix (partner just called)
#   refrac_dur  – minimum refractory duration (seconds)
#   partner_window – how long a partner's call influences transitions (seconds)
#
# Key contrast with Models 1 & 2:
#   - NO intrinsic periodicity: calls are Poisson-like within each state
#   - NO phase coupling: no sin(Δφ) or phase resets
#   - Coordination is through STATE SUPPRESSION, not timing adjustment
#   - Same individual can show very different patterns depending on context
#
# Diagnostic signatures:
#   - EXPONENTIAL / GAMMA ICI distribution (not periodic modes)
#   - Bursty, intermittent calling (clusters separated by silences)
#   - WEAK or ABSENT phase relationships at fine timescales
#   - Slow modulation over seconds-to-minutes, not cycle-by-cycle
# ============================================================================

# Packages needed to run the script:
pacman::p_load(ggplot2, dplyr, patchwork, cowplot)

# ============================================================================
# Simulation function
# ============================================================================

simulate_state_dependent <- function(
    duration = 500,              # Total simulation time (seconds)

    # ---- Calling rates per state (Hz = calls per second) ----
    # These are instantaneous rates; P(call in dt) = rate * dt.
    # Active: ~1 call/sec on average when continuously active
    # Listening: ~0.05 calls/sec (heavily suppressed)
    # Refractory: essentially 0 (but we allow a tiny rate for robustness)
    lambda_active     = 1.0,
    lambda_listening  = 0.05,
    lambda_refractory = 0.001,

    # ---- State transition probabilities (per timestep) ----
    # BASELINE transitions (when partner has NOT called recently):
    #   Active → Active:     high (tend to stay active)
    #   Active → Listening:  low  (spontaneous quieting)
    #   Active → Refractory: 0    (only enter refractory after own call)
    #   Listening → Active:  moderate (eventually return to calling)
    #   Listening → Listening: moderate
    #   Refractory → Active:  handled by refrac_dur (see below)
    p_active_to_listening_base  = 0.002,   # Spontaneous: Active → Listening
    p_listening_to_active_base  = 0.02,    # Recovery:    Listening → Active

    # PARTNER-TRIGGERED transitions (when partner called recently):
    #   The key mechanism: hearing the partner SUPPRESSES calling by
    #   pushing the focal individual from Active → Listening.
    #   This is the "probabilistic gating" that creates alternation-like
    #   patterns without any oscillator.
    p_active_to_listening_partner = 0.3,   # Suppression: Active → Listening
    p_listening_to_active_partner = 0.005, # Stay quiet longer after partner

    # ---- Refractory period ----
    # After emitting a call, the individual MUST stay in Refractory for
    # at least this many seconds. This prevents double-calling and creates
    # a minimum inter-call interval.
    refrac_dur = 0.3,

    # ---- Partner influence window ----
    # How long (seconds) after a partner's call does the partner-triggered
    # transition matrix remain active? This captures the temporal extent
    # of the suppressive effect.
    partner_window = 0.5,

    # ---- Call duration parameters ----
    d_call = 0.2,                # Mean call duration (seconds)
    duration_noise_sd = 0.1,     # SD of call duration noise

    # ---- Integration ----
    dt = 0.01                    # Timestep (seconds)
) {

  # ---- Set up time grid ----
  # Like Model 2, we use fixed-timestep simulation because state
  # transitions happen at every moment (though the model is fundamentally
  # different: no continuous phase coupling, just probabilistic state changes).
  times   <- seq(0, duration, by = dt)
  n_steps <- length(times)

  # ---- State arrays ----
  # States: 1 = Active, 2 = Listening, 3 = Refractory
  # We track states for both individuals so we can inspect dynamics later.
  states <- matrix(0L, nrow = n_steps, ncol = 2)

  # Both individuals start Active (could randomize, but Active is a
  # reasonable default -- they're about to begin an interaction).
  states[1, 1] <- 1L
  states[1, 2] <- 1L

  # ---- Output storage ----
  voc_times_1     <- numeric()
  voc_times_2     <- numeric()
  voc_durations_1 <- numeric()
  voc_durations_2 <- numeric()

  # ---- Track when each individual last called ----
  # Needed to compute: (a) partner influence window, (b) refractory exit.
  last_call_time_1 <- -Inf   # -Inf means "hasn't called yet"
  last_call_time_2 <- -Inf

  # ---- Main simulation loop ----
  for (i in 2:n_steps) {

    t_now <- times[i]

    # ================================================================
    # Process each individual (loop over ind = 1, 2)
    # ================================================================
    for (ind in 1:2) {

      current_state <- states[i-1, ind]
      partner_ind   <- 3 - ind   # if ind=1 → partner=2, and vice versa

      # When did the PARTNER last call?
      partner_last_call <- ifelse(ind == 1, last_call_time_2, last_call_time_1)
      # When did THIS individual last call? (for refractory timing)
      own_last_call     <- ifelse(ind == 1, last_call_time_1, last_call_time_2)

      # ---- Did the partner call recently? ----
      # If the partner called within the last `partner_window` seconds,
      # we use the partner-triggered transition probabilities.
      # This is the core coordination mechanism: hearing the partner
      # changes the focal individual's state dynamics.
      partner_called_recently <- (t_now - partner_last_call) < partner_window

      # ---- State transitions ----
      if (current_state == 3L) {
        # REFRACTORY state: stay here until minimum duration has passed.
        # This is a hard constraint — the animal cannot call again until
        # the refractory period is over, regardless of partner activity.
        time_in_refractory <- t_now - own_last_call

        if (time_in_refractory >= refrac_dur) {
          # Refractory period is over → transition to Active or Listening.
          # If the partner called recently, more likely to go to Listening
          # (the suppressive effect persists through the refractory period).
          if (partner_called_recently) {
            states[i, ind] <- ifelse(runif(1) < 0.6, 2L, 1L)  # 60% → Listening
          } else {
            states[i, ind] <- ifelse(runif(1) < 0.7, 1L, 2L)  # 70% → Active
          }
        } else {
          # Still in refractory — stay put
          states[i, ind] <- 3L
        }

      } else if (current_state == 1L) {
        # ACTIVE state: can transition to Listening.
        # Key mechanism: partner's call INCREASES probability of → Listening.
        if (partner_called_recently) {
          # Partner called recently → strong push toward Listening.
          # This is the "probabilistic gating" that suppresses calling.
          p_to_listen <- p_active_to_listening_partner
        } else {
          # No recent partner call → low spontaneous transition rate.
          p_to_listen <- p_active_to_listening_base
        }

        if (runif(1) < p_to_listen) {
          states[i, ind] <- 2L   # → Listening
        } else {
          states[i, ind] <- 1L   # Stay Active
        }

      } else {
        # LISTENING state: can transition back to Active.
        # Key mechanism: partner's call DECREASES probability of → Active
        # (keeps the individual quiet longer).
        if (partner_called_recently) {
          p_to_active <- p_listening_to_active_partner
        } else {
          p_to_active <- p_listening_to_active_base
        }

        if (runif(1) < p_to_active) {
          states[i, ind] <- 1L   # → Active
        } else {
          states[i, ind] <- 2L   # Stay Listening
        }
      }

      # ---- Call emission ----
      # Now that we know the current state, determine calling probability.
      # Calls are emitted stochastically (Poisson-like) with state-dependent rate.
      s <- states[i, ind]

      if (s == 1L) {
        p_call <- lambda_active * dt
      } else if (s == 2L) {
        p_call <- lambda_listening * dt
      } else {
        p_call <- lambda_refractory * dt
      }

      # Emit call?
      if (runif(1) < p_call) {
        # A call is emitted!
        call_dur <- max(0.01, d_call + rnorm(1, 0, duration_noise_sd))

        if (ind == 1) {
          voc_times_1     <- c(voc_times_1, t_now)
          voc_durations_1 <- c(voc_durations_1, call_dur)
          last_call_time_1 <- t_now
        } else {
          voc_times_2     <- c(voc_times_2, t_now)
          voc_durations_2 <- c(voc_durations_2, call_dur)
          last_call_time_2 <- t_now
        }

        # After calling → enter Refractory (mandatory).
        # This is deterministic: every call triggers a cooldown period.
        states[i, ind] <- 3L
      }
    }
  }

  # ---- Return results ----
  list(
    voc_times_1     = voc_times_1,
    voc_times_2     = voc_times_2,
    voc_durations_1 = voc_durations_1,
    voc_durations_2 = voc_durations_2,
    # Full state trajectories for inspection
    times  = times,
    states = states,    # n_steps x 2 matrix (1=Active, 2=Listening, 3=Refractory)
    params = list(
      lambda_active = lambda_active,
      lambda_listening = lambda_listening,
      lambda_refractory = lambda_refractory,
      p_active_to_listening_base = p_active_to_listening_base,
      p_listening_to_active_base = p_listening_to_active_base,
      p_active_to_listening_partner = p_active_to_listening_partner,
      p_listening_to_active_partner = p_listening_to_active_partner,
      refrac_dur = refrac_dur,
      partner_window = partner_window,
      d_call = d_call,
      duration_noise_sd = duration_noise_sd,
      dt = dt,
      duration = duration
    )
  )
}
