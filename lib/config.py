"""
    config module
    ~~~~~~~~~~~~~

    Implements the configuration related objects.
"""

import types
import errno
from collections import UserDict


class Config(UserDict):
    """Instantiate a configuration object in several ways.

    :param default: an optional dictionary of default values."""

    def __init__(self, default=None):
        self.data = {}
        if default is not None:
            try:
                self.update(default)
            except TypeError:
                pass

    def from_pyfile(self, filename, silent=False):
        """Updates the values in config from a Python file. This function
        behaves as if the file was imported as a module with the
        :meth `from_object` function.
        :param filename: the filename of the config.
        :param silent: set to ``True`` if you want silent failure for missing
                       files.
        """
        filename = filename
        d = types.ModuleType('config')
        d.__file__ = filename
        try:
            with open(filename, mode='rb') as config_file:
                exec(compile(config_file.read(), filename, 'exec'), d.__dict__)
        except IOError as e:
            if silent and e.errno in (errno.ENOENT, errno.EISDIR):
                return False
            e.strerror = 'Unable to load configuration file (%s)' % e.strerror
            raise
        self.from_object(d)
        return True

    def from_object(self, obj):
        """Updates the values from the given object.  An object can be of one
           of the following two types:

           -  a string: in this case the object with that name will be imported
           -  an actual object reference: that object is used directly

        Objects are usually either modules or classes. :meth:`from_object`
        loads only the uppercase attributes of the module/class. A ``dict``
        object will not work with :meth:`from_object` because the keys of a
        ``dict`` are not attributes of the ``dict`` class.

        Example of module-based configuration::
            app.config.from_object('yourapplication.default_config')
            from yourapplication import default_config
            app.config.from_object(default_config)

        You should not use this function to load the actual configuration but
        rather configuration defaults.  The actual config should be loaded
        with :meth:`from_pyfile`.

       :param obj: an import name or object
       """

        for key in dir(obj):
            if key.isupper():
                self.data[key] = getattr(obj, key)

    def __setitem__(self, key, value):
        self.data[str(key)] = value




