import streamlit as st
from dolfin import *
from PIL import Image
 
menu = st.sidebar.radio("***",
    ("Библиотеки Python", 
    "Текст программы", 
    "Результат работы",
    "Параметрические расчеты",    
    "streamlit - Parameters",
    "streamlit - Plot",   
    )
)
  
if menu == "Библиотеки Python":
    r"""
##### Библиотеки Python
* fenics - вычислительная платформа конечно-элементного анализа
* matplotlib - графика

    """
      
if menu == "Текст программы":
    r"""
##### Текст программы

**pr.py**

    """
    code = """
from dolfin import *
import matplotlib.pyplot as plt

# Parameters
n = 20
eps = Constant(0.001)

# Domain and mesh
mesh = UnitSquareMesh(n, n)
plot(mesh, title="mesh")
plt.savefig("gr1.png", format="png", dpi=200)
plt.clf()

# Define function space
V = FunctionSpace(mesh, "CG", 2)

# Subdomain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] + x[1] > 1.0 -  DOLFIN_EPS and on_boundary
# Define boundary condition
bc = DirichletBC(V, Constant(0.0), DirichletBoundary())

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("x[0]", degree=2)
a = eps*inner(grad(u), grad(v))*dx + u*v*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)
print ("max|u(x)| = ", u.vector().norm('linf') )

# Plot  
plot(u, title="Solution")
plt.savefig("gr2.png", format="png", dpi=200)
plt.clf()
"""
    st.code(code, language="python")  
          
if menu == "Результат работы":
    r"""
##### Результат работы

Параметры
* $n = 20$
* $\varepsilon = 0.001$

Консольный вывод
    """
    code = """
Solving linear variational problem.
max|u(x)| =  0.8571756975999396
    """
    st.code(code, language="bash")  
    
    r"""
Графики: gr1.png, gr2.png
    """
    c1, c2 = st.columns(2)
    image1 = Image.open("pages/figs/gr1.png")
    c1.image(image1)
    image2 = Image.open("pages/figs/gr2.png")
    c2.image(image2)
    
if menu == "Параметрические расчеты":

    import matplotlib.pyplot as plt
    import numpy as np

    r"""
##### Параметрические расчеты
* Равномерная сетка $n \times n$    
* Малый параметр $\varepsilon = 10^{-p}$
    """
        
    # Parameters    
    c1, cc, c2 = st.columns([5,1,5])
    n = c1.selectbox("N", [10, 20, 40], index = 1)
    pa = range(1,10)
    p = c2.select_slider("p", options=pa, value=3)
    eps = Constant(10**(-p))
    
    # Domain and mesh
    mesh = UnitSquareMesh(n, n)
    
    # Define function space
    V = FunctionSpace(mesh, "CG", 2)
    
    # Subdomain for Dirichlet boundary condition
    class DirichletBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] + x[1] > 1.0 -  DOLFIN_EPS and on_boundary
    # Define boundary condition
    bc = DirichletBC(V, Constant(0.0), DirichletBoundary())
    
    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Expression("x[0]", degree=2)
    a = eps*inner(grad(u), grad(v))*dx + u*v*dx
    L = f*v*dx
    
    # Compute solution
    u = Function(V)
    solve(a == L, u, bc)
   
    # Plot  
    plot(u, title="Solution")
    st.pyplot(plt.gcf())
    plt.clf()
    str = "max|u(x)| = " + str(u.vector().norm('linf'))
    st.write(str)


if menu == "streamlit - Parameters":
    r"""
##### streamlit - Parameters
    """

    code = """
    # Parameters    
    c1, cc, c2 = st.columns([5,1,5])
    n = c1.selectbox("N", [10, 20, 40], index = 1)
    pa = range(1,10)
    p = c2.select_slider("p", options=pa, value=3)
    eps = Constant(10**(-p))
    """
    st.code(code, language="python")  
    
if menu == "streamlit - Plot":
    r"""
##### streamlit - Plot
    """

    code = """
    # Parameters    
    # Plot  
    plot(u, title="Solution")
    st.pyplot(plt.gcf())
    plt.clf()
    str = "max|u(x)| = " + str(u.vector().norm('linf'))
    st.write(str)
    """
    st.code(code, language="python")  


    
    
       

