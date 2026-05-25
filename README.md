# CalGNSSC

Pure C GNSS LLA to J2000 state fitting and short-term prediction project.

## Build

Open `CalGNSSC.sln` in Visual Studio and build `Release|x64`. The solution now
contains a standalone `calgnss_orbit_mini` static library project plus a `CalGNSSC`
console test entry that references the library. Or use CMake on other
platforms:

```sh
cmake -S . -B build
cmake --build build --config Release
```

The core code uses only the C standard library and `math.h`. The Visual Studio
project is only a build wrapper; the same `include/calgnss_orbit_mini.h` and
`src/calgnss_orbit_mini.c` can be reused in Linux or RTOS builds.

## C API

`cg_context_t` is opaque. Create it through `cg_context_create`, feed ordered
LLA observations with `cg_context_push`, query J2000 position/velocity with
`cg_context_query_state`, and release it with `cg_context_destroy`.

```c
cg_observation_t buffer[CG_DEFAULT_OBSERVATION_CAPACITY];
cg_context_t *context = NULL;
cg_state_t state;

status = cg_context_create(&context, buffer, CG_DEFAULT_OBSERVATION_CAPACITY, &options);
status = cg_context_push(context, &observation);
status = cg_context_query_state(context, &query_time, &state);
cg_context_destroy(context);
```

## Default Test Entry

From this directory:

```sh
x64\Release\CalGNSSC.exe
```

The test entry reads its export settings directly from
`examples/export_j2000.c`:

```text
k_input_csv = "..\20260518\67113_LLA_Position.csv"
k_output_csv = "..\20260518\67113_J2000_Calculated.csv"
k_start_time_text = "18 May 2026 04:00:00.000"
k_end_time_text = "19 May 2026 04:00:00.000"
k_step_seconds = 1.0
```

## Method

- WGS-84 geodetic LLA is converted to ECEF.
- ECEF is transformed to J2000/GCRS with a compact IAU 2000B
  precession-nutation and Greenwich apparent sidereal time model.
- In-span queries use a local Chebyshev least-squares position fit; velocity is
  the analytical derivative of that fit.
- Future queries use RK4 J2 propagation from the final state.
  When at least one prior orbit of history exists, the propagated position is
  corrected with the previous-orbit residual at the same horizon, which fixes
  the unstable short extrapolation behavior seen in the Python prototype.
