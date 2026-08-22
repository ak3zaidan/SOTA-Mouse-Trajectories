# pointer_motion.py

> This is the highest quality cursor generator you can use

Generates human-like mouse trajectories between a series of waypoints. Instead of moving in straight lines at constant speed, the generated paths curve, accelerate, decelerate, and jitter the way a real hand does.

## Usage

```python
from pointer_motion import traj

# coordinates = [[x0, x1, ...], [y0, y1, ...]]  — waypoints to travel through
move, move_coalesced = traj(
    coordinates,      # waypoint x/y lists
    screen_w,         # screen width  (px)
    screen_h,         # screen height (px)
    160,              # sampling frequency (Hz)
    1,                # sv — noise scale factor
    True,             # coal — also return coalesced (60 Hz batched) events
    True,             # remove_v0 — drop stationary samples
)
```

### Output

`move` is a list of 5 parallel arrays, one entry per stroke (waypoint pair):

| Index | Contents |
|-------|----------|
| `move[0]` | x positions (px, rounded) |
| `move[1]` | y positions (px, rounded) |
| `move[2]` | timestamps (ms, relative to start) |
| `move[3]` | ideal velocity at each point |
| `move[4]` | observed velocity (recomputed from final noisy pixels) |

With `coal=True` a second structure is returned where the samples are grouped into batches the way an OS dispatches pointer events (~every 16.7 ms at 60 Hz).

## How it works

Each pair of consecutive waypoints becomes one **stroke**, built in these stages:

1. **Log-normal velocity profile.** Human pointer strokes follow a log-normal speed curve: quick ramp-up to a peak, then a long tail of deceleration. Each stroke gets a peak velocity and duration scaled to its distance, and an iterative solver (`errorf`) tunes the log-normal parameters (`mu`, `sigma`) until the stroke lands on target with a near-zero end velocity.

2. **Curved approach angle.** The stroke doesn't head straight at the target. It starts at a randomized angle offset from the direct line (`ang_s`), and the direction is swept toward the end angle over time using an error-function blend — producing a natural arc rather than a straight shot.

3. **End alignment.** Past the velocity peak, the tail of the log-normal curve is replaced with a straight-line velocity decay aimed directly at the target, so the pointer settles exactly on it. Any remaining pixel error is distributed proportionally across the stroke's points.

4. **Bezier smoothing.** A handful of random support points are picked along the path and quadratic Bezier segments are fit between them, each bowed slightly to a random side. This adds the gentle, alternating wobble of real hand movement.

5. **Noise.** Small velocity-dependent jitter is added, split into tangential (along the movement) and perpendicular components — faster movement gets more noise, and the start/end of each stroke stay clean.

6. **Post-processing.** Positions are clamped to the screen, rounded to integer pixels, and timestamps converted to ms. Optionally, zero-movement samples are removed (`remove_v0`) and a random 0.5–1.1 s idle pause is appended after each stroke to simulate the hesitation before a click.

## Helper functions

- `biased_random(mean, std, lo, hi)` — samples a normal distribution, rejecting values outside bounds. Used for the start-angle randomization.
- `add_coalesced_events(move)` — groups raw samples into 60 Hz event batches, mimicking OS pointer-event coalescing.
- `filter_v0(move)` — removes samples where the pointer didn't move.
- `add_click_delay(move, i, frequ)` — appends stationary points after a stroke to model the pause before clicking.
