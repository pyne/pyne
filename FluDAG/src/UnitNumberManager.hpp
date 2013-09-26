class UnitNumberManager
{
  public:
    /**
     * \brief Constructor
     */
    UnitNumberManager();

    static int BAD_UNIT_NUMBER = 0;

    int getUnitNumber(std::string name);

  private:
    /* Start and end logical (Fortran-style) unit numbers for record writing */
    static int START_UNIT = -21;
    static int END_UNIT   = -99;
    /* The number of types of cards to be written. Card 'types' are distinguished by
     * the combination of particle and tally type 
     */
    static int num_units_in_use = 0;

    std::map<std::string, int> nameNumberMap;
    int getNextUnitNumber();
}
