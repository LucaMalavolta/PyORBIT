def pytrades():
    pass

def DiffEvol():
    pass

def PyPolyChord():
    pass

def PolyChord():
    pass

def celerite():
    pass

def ttvfast():
    pass

def george():
    pass

def batman():
    pass

## absurd workaround to fix the lack of celerite in the system
def Celerite_QuasiPeriodicActivity():
    pass


class dummy_one:
    def __init__(self):
        self.terms = dummy_two()

class dummy_two:
    def __init__(self):
        self.Term = dummy_three(0, 0, 0, 0)
    #def Term(self):
    #    return

class dummy_three:
    def __init__(self, a, b, c, d):
        self.Term = 0
