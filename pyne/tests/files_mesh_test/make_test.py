from pyne.mesh import Mesh

my_mesh = Mesh(structured_coords=[[-1,0,1],[-1,0,1],[0,1]], structured=True)
volumes1 = list(my_mesh.structured_iterate_hex("xyz"))
volumes2 = list(my_mesh.structured_iterate_hex("xyz"))
flux_tag = my_mesh.mesh.createTag("flux", 1, float)
error_tag = my_mesh.mesh.createTag("flux_error", 1, float)
flux_data = [1.1, 2.2, 3.3, 4.4]
error_data = [0.1, 0.2, 0.3, 0.4]
flux_tag[volumes1] = flux_data
error_tag[volumes2] = error_data
my_mesh.mesh.save("test_mesh_2x2x1_2.h5m")
