import collections

"""
A custom dictionary object with a recursive update.
"""

class CustomDict(dict):
    def update(self, u):
        for k, v in u.iteritems():
            if isinstance(self, collections.Mapping):
                if isinstance(v, collections.Mapping):
                    r = self.get(k, CustomDict({})).update(v)
                    self[k] = r
                else:
                    self[k] = u[k]
            else:
                self = {k: u[k]}
        return self
