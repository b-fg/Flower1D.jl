# Abstract element
abstract type AbstractElement{T <: AbstractFloat} end

"""
    StandardElement{T <: AbstractFloat} <: AbstractElement{T}

Standard (aka reference) element structure. It stores properties and variables of standard elements.

# Fields
- `type :: Integer`. Type of standard element. Currently only 1-D element (line) supported.
- `order :: Integer`. Order of the polynomials. Currently only implemented for p = 1 to 5.
- `quadratureType :: String`. Type of quadrature used in the element. Supports `"equidistant"` and `"legendre"`.
- `innerPoints :: Vector{T}`. 1-D array storing the location of the inner (aka integration) points of the element. Size = p + 1.
- `facePoints :: Vector{T}`. 1-D array storing the location of the two face integration points of the element. Size = 2.
- `quadratureWeights :: Vector{T}`. 1-D array storing weights of the Gaussian quadrature at inner points. Size = 2.
- `basisInner :: Matrix{T}`. 2-D array storing the values of the Lagrange polynomials at each inner point.
    (a bit pointless to store it, just for debugging). Size = p + 1, p + 1
- `basisFace :: Matrix{T}`. 2-D array storing the values of the Lagrange polynomials at each face. Size = 2, p + 1.

# Arguments
- `type :: Integer`. (As above).
- `order :: Integer`. (As above).
- `quadratureType :: String`. (As above).
- `T :: Type`. Float precision type.
"""
struct StandardElement{T <: AbstractFloat} <: AbstractElement{T}
    type :: Integer
    order :: Integer
    quadratureType :: String
    innerPoints :: Vector{T}
    facePoints :: Vector{T}
    quadratureWeights :: Vector{T}
    basisInner :: Matrix{T}
    basisFace :: Matrix{T}

    function StandardElement(type :: Integer, order :: Integer, quadratureType :: String; T = Float64)
        innerPoints, weights = quadPointsAndWeights(type, order, quadratureType)
        facePoints = [-1.0, 1.0]
        basisInner = zeros(T, length(innerPoints), length(innerPoints))
        basisFace = zeros(T, length(facePoints), length(innerPoints))
        for i ∈ 1:length(innerPoints)
            basisInner[i, :] = lagrangeBasisAt(innerPoints[i], innerPoints)
        end
        for i ∈ 1:length(facePoints)
            basisFace[i, :] = lagrangeBasisAt(facePoints[i], innerPoints)
        end
        new{T}(type, order, quadratureType, innerPoints, facePoints, weights, basisInner, basisFace)
    end
end

"""
    Element{T <: AbstractFloat} <: AbstractElement{T}

Element structure. It stores properties and variables of elements in the physical space.

# Fields
- `vertexPoints :: Vector{T}`. 1-D array storing the geometrical nodes that define the element.
- `innerPoints :: Vector{T}`. 1-D array storing the physical location of the inner (aka integration) points.
- `facePoints :: Vector{T}`. 1-D array storing the physical location of the face integration points.
    Note that for 1-D it is equivalent to `vertexPoints`.
- `Δx :: AbstractFloat`. Float value storing the size of the element.

# Arguments
- `vertexPoints :: Vector{T}`. (As above).
- `standardElement :: StandardElement`. Standard element struct..
- `T :: Type`. Float precision type.
"""
struct Element{T <: AbstractFloat} <: AbstractElement{T}
    vertexPoints :: Vector{T} # geometrical mesh nodes defining the element
    innerPoints :: Vector{T} # solution points within the element (according to quadrature)
    facePoints :: Vector{T} # flux points located at the faces
    Δx :: AbstractFloat # element size

    function Element(vertexPoints :: Vector{TF} where TF<: AbstractFloat, standardElement :: StandardElement; T = Float64)
        innerPoints = xFromξ.(standardElement.innerPoints, vertexPoints[1], vertexPoints[end])
        facePoints = xFromξ.(standardElement.facePoints, vertexPoints[1], vertexPoints[end])
        Δx = vertexPoints[end] - vertexPoints[1]
        new{T}(vertexPoints, innerPoints, facePoints, Δx)
    end
end

"""
    PolynomialData{T <: AbstractFloat}

Polynomial data required to compute RHS based on standard element data.

# Fields
- `Dij :: Matrix{T}`. 2-D array storing the Lagrange derivative matrix: first-order derivative of polynomials at inner points.
- `correctionDerivativeL :: Vector{T}`. 1-D array storing the first-order derivative of the left correction function.
- `correctionDerivativeR :: Vector{T}`. 1-D array storing the first-order derivative of the right correction function.

# Arguments
- `innerPoints :: Vector{T}`. 1-D array storing the inner points location in the standard element.
- `P :: Integer`. Element polynomial order.
- `C :: Union{Integer, AbstractFloat}`. Integer or Float value for the free parameter of the vcjh correction function.
    * C = 0::Int ➡ DG nodal
    * C = 1::Int ➡ Spectral Difference (SD)
    * C = 2::Int ➡ Huynh
    * C <: AbstractFloat, this exact value is used.
"""
struct PolynomialData{T <: AbstractFloat}
    Dij :: Matrix{T}
    correctionDerivativeL :: Vector{T}
    correctionDerivativeR :: Vector{T}

    function PolynomialData(innerPoints :: Vector{TF} where TF<: AbstractFloat, P :: Integer, C :: Union{Integer, AbstractFloat};
        T = Float64)
        Dij = lagrangeDerivativeMatrix(innerPoints)
        correctionDerivative = vcjh(P, innerPoints; C = C, derivativeChoice = true)
        new{T}(Dij, correctionDerivative[:, 1], correctionDerivative[:, 2])
    end
end

"""
    Mesh{T <: AbstractFloat}

Mesh structure storing elements and polynomial to compute the RHS of the PDEs, among other variables.

The mesh indexing is structured as follows (example with 3 elements and 2 ghost elements, total of 5 elements):

faces:      na    1     2     3     4     na
            x--G--x-----x-----x-----x--G--x
elements:      1     2     3     4     5

# Fields
- `standardElement :: StandardElement`. Standard element of the mesh. Since hybrid meshes are not supported, a single instance is stored.
- `polynomialData :: PolynomialData`. Polynomial data required to compute the RHS.
- `elements :: Vector{Element}`. 1-D array storing the `Element`s of the mesh.
- `elementsNodes :: Matrix{T}`. 2-D array storing geometrical nodes of the mesh that define each element [element index, node index].
- `JijInv :: Vector{T}`. 1-D array the inverse of each element Jacobian.
- `elementsType :: Integer`. Integer defining the type of element. Only `1` (line) supported.
- `elementsOrder :: Integer`. Integer defining the order of the elements.
- `C :: Union{Integer, AbstractFloat}`. Integer or Float value for the free parameter of the vcjh correction function. See `PolynomialData` for details.
- `startIndex :: Integer`. Index of the first non-ghost element.
- `endIndex :: Integer`. Index of the last non-ghost element.

# Arguments
- `elementsType :: Integer`. See above.
- `elementsOrder :: Integer`. See above.
- `quadratureType :: String`. Type of quadrature used in the element. Supports `"equidistant"` and `"legendre"`.
- `elementsNodes :: Matrix{T}`. See above.
- `C :: Union{Integer, AbstractFloat}`. See above.
- `T :: Type`. Float precision type.
"""
struct Mesh{T <: AbstractFloat}
    standardElement :: StandardElement
    polynomialData :: PolynomialData
    elements :: Vector{Element}
    elementsNodes :: Matrix{T}
    JijInv :: Vector{T}
    elementsType :: Integer
    elementsOrder :: Integer
    C :: Union{Integer, AbstractFloat}
    startIndex :: Integer
    endIndex :: Integer

    function Mesh(elementsType :: Integer, elementsOrder :: Integer, quadratureType :: String,
        elementsNodes :: Matrix{TF} where {TF}, C :: Union{Integer, AbstractFloat}; T = Float64)
        N = size(elementsNodes, 1) # number of elements

        standardElement = StandardElement(elementsType, elementsOrder, quadratureType; T = T)
        polynomialData = PolynomialData(standardElement.innerPoints, elementsOrder, C; T = T)

        elements = Vector{Element}(undef, N + 2)
        JijInv = Vector{T}(undef, N)
        for n ∈ 1:N
            elements[n + 1] = Element(elementsNodes[n, :], standardElement; T = T)
            JijInv[n] = 2.0 / (elementsNodes[n, 2] - elementsNodes[n, 1])
        end
        # Add ghost elements at start and end of elements Vector
        elements[1] = elements[N + 1]
        elements[end] = elements[2]

        new{T}(standardElement, polynomialData, elements, elementsNodes, JijInv, elementsType, elementsOrder, C, 2, N + 1)
    end
end

meshSize(mesh :: Mesh) = length(mesh.elements)
meshSizeNoGhosts(mesh :: Mesh) = length(mesh.elements) - 2
meshOrder(mesh :: Mesh) = mesh.elementsOrder
meshDOF(mesh :: Mesh) = meshSize(mesh) * (mesh.elementsOrder + 1)
meshDOFNoGhost(mesh :: Mesh) = meshSizeNoGhosts(mesh) * (mesh.elementsOrder + 1)
meshΔxMin(mesh :: Mesh) = minimum([mesh.elements[n].Δx for n ∈ 1:meshSize(mesh)])
meshΔx(mesh :: Mesh) = [mesh.elements[n].Δx for n ∈ 2:meshSizeNoGhosts(mesh) + 1]
element2faces(elementIndex :: Integer) = (elementIndex - 1, elementIndex)
face2elements(faceIndex :: Integer) = (faceIndex, faceIndex + 1)
startIndex(mesh :: Mesh) = mesh.startIndex
endIndex(mesh :: Mesh) = mesh.endIndex

"""
    integrationPointsLocation(mesh :: Mesh)

Returns a `Vector` containing the location of the integration points of all elements in the mesh in a flattened style.
"""
function integrationPointsLocation(mesh :: Mesh)
    x = Vector{eltype(mesh.elementsNodes)}(undef, meshDOFNoGhost(mesh))
    k = meshOrder(mesh) + 1 # number of solution points per element
    for i in 2:meshSize(mesh) - 1
        x[(i - 2) * k + 1:(i - 2) * k + k] = mesh.elements[i].innerPoints[:]
    end
    return x
end