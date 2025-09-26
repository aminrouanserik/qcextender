import sxs

print(sxs.sxs_directory("config"))  # Probably returns some path in your home directory
sxs.write_config(download=True, cache=True, auto_supersede=False)
print(sxs.read_config())
