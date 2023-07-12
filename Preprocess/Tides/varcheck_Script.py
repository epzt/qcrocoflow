def varcheck(namespace, varname, default_value):
    # Check if the variable exists in the namespace
    var = namespace.get(varname, default_value)

    # If the variable is empty, set it to the default value
    if var is None:
        var = default_value

    # Set the variable in the namespace
    namespace[varname] = var


env = {}
varcheck(env, 'my_var', 'default_value')
print(env['my_var'])
