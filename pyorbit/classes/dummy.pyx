def pytrades():
    pass

def pyde():
    pass

def PolyChord():
    pass

def celerite():
    pass

def ttvfast():
    pass

def george():
    pass

## absurd workaround to fix the lack of celerite in the system

class dummy_one:
    def __init__(self):
        self.terms = dummy_two()

class dummy_two:
    def __init__(self):
        self.Term = dummy_three(0, 0, 0)
    #def Term(self):
    #    return

class dummy_three:
    def __init__(self, a, b, c):
        self.Term = 0
