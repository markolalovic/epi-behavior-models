# Local SQL databases for ABC-SMC logging

- SQL databases `sqlite` files used to log and store intermediate results from ABC-SMC parameter inference runs. 

- These are used by pyABC framework as a storage system, see [pyABC data store documentation](https://pyabc.readthedocs.io/en/latest/datastore.html).

- These files are not tracked by version control, see `.gitignore`, as they can be large (~60 MB per location, and ~300 MB for model selection results).

## Example
To configure pyABC to log ABC-SMC runs into "test.db" stored in this directory, specify the path to the local database:

```python
import os
import pyabc

local_db_path = "../data/local-db/"
db_path = os.path.join(local_db_path, "test.db")
```

And initialize pyABC logging:
```python
abc = pyabc.ABCSMC(...)
abc.new("sqlite:///" + db_path, {"data": obs_deaths})
```

