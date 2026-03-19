class NeurobridgeChecks {

    static void requireParam(value, String name) {
        if( value == null || value.toString().trim() == '' ) {
            throw new IllegalArgumentException("Missing --${name}")
        }
    }
}