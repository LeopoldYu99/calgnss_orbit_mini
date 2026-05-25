# calgnss_orbit_mini

Standalone C library for GNSS LLA to J2000 state fitting and short-term
prediction.

## Build

Use CMake from this directory:

```sh
cmake -S . -B build
cmake --build build --config Release
```

The project builds both a shared library for Python `ctypes` loading and a
static library for native C/C++ integration:

```text
build/Release/calgnss_orbit_mini.dll
build/Release/calgnss_orbit_mini.lib
build/Release/calgnss_orbit_mini_static.lib
```

The core code uses only the C standard library and `math.h`; on Unix-like
systems CMake links `libm` automatically.

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

The in-span fit window is fixed inside the library at 60 minutes. To limit how
many recent observations can participate in a query, choose the size of the
buffer passed to `cg_context_create`.

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
