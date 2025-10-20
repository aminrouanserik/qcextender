# Not a test file, rather a set-up file (change?)

import sxs

print(sxs.sxs_directory("config"))  # Probably returns some path in your home directory
sxs.write_config(
    download=True, cache=True, auto_supersede=True, ignore_deprecation=True
)
print(sxs.read_config())
