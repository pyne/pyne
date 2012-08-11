



rxr = simulation.definition.SystemDefinition()

rxr.add_material()

rxr.add_material()

rxr.add_surface()

rxr.add_surface()

rxr.add_cell()

rxr.add_cell()

opts = simulationdef.OptionsDefinition()

crit_source_name = opts.add_criticality_source()

opts.add_criticality_points()

opts.add_cellflux_tally()

rxr.save("system1")

opts.save("options1")

enrichments = [0.01 0.02 0.03]

for this_enrich in enrichments:
    rxr.material['UO2'] = material(blah blah)
    inp = MCNPInput("input1", rxr, opts)
    inp.write()

    inp.system.material['UO2'] = material(blah blah)
    inp.options.card[crit_source_name].cycles = 1000
    inp.write()
    



