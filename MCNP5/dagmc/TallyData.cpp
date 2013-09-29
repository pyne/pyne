
double& MeshTally::get_data(std::vector<double>& data,
                            moab::EntityHandle tally_point,
                            unsigned energy_bin)
{
    assert(energy_bin < num_energy_bins);
    int index = get_entity_index(tally_point) * num_energy_bins + energy_bin;
    return data[index];
}

