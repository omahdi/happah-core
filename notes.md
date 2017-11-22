## Notes

### Converting OFF to HPH

```
#include "happah/format.h"

auto content = format::off::read("prism.off");
auto mesh = make_triangle_mesh<VertexP3>(content);
auto graph = make_triangle_graph(mesh);
auto surface0 = make_spline_surface(graph);
auto surface1 = elevate(surface0);

write(mesh, "prism.tm.p3.s.hph");
write(surface1, "prism.ss5.bz.3.hph");
```

### Creating a Sunny Projective Structure

```
auto path = undegenerate(graph, trim(graph, cut(graph)));
validate_cut(graph, path);
auto analysis = analyze(graph, path);
auto structure = make_projective_structure(std::get<0>(analysis), std::get<2>(analysis));

auto border = Indices();
border.reserve(std::get<0>(analysis).size());
for(auto i = hpindex(1), end = hpindex(1 + 3 * std::get<0>(analysis).size()); i != end; i += 3) border.push_back(i);
auto sun = detail::make_sun(std::get<0>(analysis));
auto mesh = make_triangle_mesh(structure, border, 0, Point3D(0, 0, 1), Point3D(std::get<0>(sun)[0], 1), Point3D(std::get<0>(sun)[1], 1));
```
