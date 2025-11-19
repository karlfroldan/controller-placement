from dotenv import dotenv_values

def ampl_license() -> str:
    env = dotenv_values('.env')

    license = env['AMPL_LICENSE_UUID']
    if not license or license == '':
        raise Exception('`AMPL_LICENSE_UUID` is not defined')
    
    return env['AMPL_LICENSE_UUID']

