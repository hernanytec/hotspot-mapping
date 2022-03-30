import abc

class BaseEstimator(metaclass=abc.ABCMeta):

    @abc.abstractclassmethod
    def fit(self):
        return

    @abc.abstractclassmethod
    def score(self):
        return