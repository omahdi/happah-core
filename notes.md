## Notes

### Converting OFF to HPH

```
#include "happah/format.h"

auto content = format::off::read("prism.off");
auto mesh0 = make_triangle_mesh(content);
auto mesh1 = TriangleMesh<VertexP3, Format::DIRECTED_EDGE>(mesh0);
auto surface0 = make_spline_surface(mesh1);
auto surface1 = elevate(surface0);

write(mesh0, "prism.tm.p3.s.hph");
write(surface1, "prism.ss5.bz.3.hph");
```

