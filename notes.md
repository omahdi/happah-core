## Notes

### Converting OFF to HPH

```
#include "happah/writers.h"

auto output = ReaderOFF::read("prism.off");
auto mesh0 = ReaderOFF::toTriangleMesh(output);
write(mesh0, "prism.tm.3.simple.hph");
auto mesh1 = TriangleMesh<VertexP3, Format::DIRECTED_EDGE>(mesh0);
auto surface0 = make_quartic_surface_spline_bez(mesh1);
auto surface1 = elevate(surface0);
write(surface1, "prism.ss5.bz.3.hph");
```

