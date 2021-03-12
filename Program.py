import NXOpen.GeometricUtilities
import NXOpen.Features
import NXOpen.GeometricAnalysis
import NXOpen
import NXOpen.Assemblies
import NXOpen.UF
import time


####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
# Settings
COLOR = 187
REFERENCE_SET = "DETAIL-NUMBERS"
LAYER = 1
FONT = "Modern"
LENGTH = 0.900
HEIGHT = 0.400
SCRIPT = NXOpen.Features.TextBuilder.ScriptOptions.Oem
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################

THE_SESSION = NXOpen.Session.GetSession()

THE_UF_SESSION = NXOpen.UF.UFSession.GetUFSession()

THE_UISESSION = NXOpen.UI.GetUI()

UF_SOLID_TYPE = NXOpen.UF.UFConstants.UF_solid_type

UF_UI_SEL_FEATURE_BODY = NXOpen.UF.UFConstants.UF_UI_SEL_FEATURE_BODY

UF_CSYS_ROOT_COORDS = NXOpen.UF.UFConstants.UF_CSYS_ROOT_COORDS

UF_CSYS_ROOT_WCS_COORDS = NXOpen.UF.UFConstants.UF_CSYS_ROOT_WCS_COORDS


def display_part(set_display=None) -> NXOpen.Part:
    if set_display is None:
        return THE_SESSION.Parts.Display
    THE_SESSION.Parts.SetDisplay(
        set_display, False, False)
    return None


def work_part(set_work=None) -> NXOpen.Part:
    if set_work is None:
        return THE_SESSION.Parts.Work
    THE_SESSION.Parts.SetWork(set_work)
    return None


def write_line(message) -> None:
    lw = NXOpen.Session.GetSession().ListingWindow
    lw.Open()
    if isinstance(message, int) or isinstance(message, float) or isinstance(message, type):
        lw.WriteLine(str(message))
        return
    if isinstance(message, NXOpen.Point3d) or isinstance(message, NXOpen.Vector3d):
        lw.WriteLine(str([message.X, message.Y, message.Z]))
        return
    if isinstance(message, list):
        lw.WriteLine(str(message))
        return
    lw.WriteLine(message)


def select_solid_bodies():
    mask = NXOpen.SelectionMaskTriple_Struct(
        UF_SOLID_TYPE, 0, UF_UI_SEL_FEATURE_BODY)
    message = "Please select tool bodies"
    title = "Selection"
    scope = NXOpen.SelectionSelectionScope.AnyInAssembly
    action = NXOpen.SelectionSelectionAction.ClearAndEnableSpecific
    include_features = False
    keep_highlighted = False
    mask_array = [mask]
    return THE_UISESSION.SelectionManager.SelectTaggedObjects(
        message, title, scope, action, include_features, keep_highlighted, mask_array)


def select_solid_body():
    mask = NXOpen.SelectionMaskTriple_Struct(
        UF_SOLID_TYPE, 0, UF_UI_SEL_FEATURE_BODY)
    message = "Please select tool bodies"
    title = "Selection"
    scope = NXOpen.SelectionSelectionScope.AnyInAssembly
    action = NXOpen.SelectionSelectionAction.ClearAndEnableSpecific
    include_features = False
    keep_highlighted = False
    mask_array = [mask]
    return THE_UISESSION.SelectionManager.SelectTaggedObject(
        message, title, scope, action, include_features, keep_highlighted, mask_array)


def get_detail_name(component: NXOpen.Assemblies.Component) -> str:
    display_name = component.DisplayName
    length = len(display_name)
    return display_name[length - 3: length]


def nxopen_point3d_to_list(point: NXOpen.Point3d) -> list:
    return [point.X, point.Y, point.Z]


def nxopen_vector3d_to_list(vector: NXOpen.Vector3d) -> list:
    return [vector.X, vector.Y, vector.Z]


def compute_distance(face1, face2, tolerance: float = .001) -> float:
    """
    Calculates the distance between to planar faces.

    Parameters
    -----
    face1 : NXOpen.Face
        The first face
    face2 : NXOpen.Face
        The second face

    Optional Parameters
    -----
    tolerance : float
        The tolerance, default value: .001
    """
    pnt1 = face1.GetEdges()[0].GetVertices()[0]
    pnt_on_plane = face2.GetEdges()[0].GetVertices()[0]
    plane_normal = normal_vector(face1)
    return THE_UF_SESSION.Vec3.DistanceToPlane(
        nxopen_point3d_to_list(pnt1),
        nxopen_point3d_to_list(pnt_on_plane),
        nxopen_vector3d_to_list(plane_normal),
        tolerance)


def start_point_edge(edge: NXOpen.Edge) -> NXOpen.Point3d:
    return edge.GetVertices()[0]


def end_point_edge(edge: NXOpen.Edge) -> NXOpen.Point3d:
    return edge.GetVertices()[1]


def normal_vector(face: NXOpen.Face) -> NXOpen.Vector3d:
    tup = THE_UF_SESSION.Modeling.AskFaceData(face.Tag)
    return NXOpen.Vector3d(tup[2][0], tup[2][1], tup[2][2])


def unitize_vector(vector: NXOpen.Vector3d, tolerance: float = .001) -> NXOpen.Vector3d:
    try:
        tup = THE_UF_SESSION.Vec3.Unitize(
            [vector.X, vector.Y, vector.Z], tolerance)
        return NXOpen.Vector3d(tup[1][0], tup[1][1], tup[1][2])
    except NXOpen.NXException:
        write_line("Error: " + str([vector.X, vector.Y, vector.Z]))


def is_parallel(vector1, vector2, tolerance=.001) -> bool:
    result = THE_UF_SESSION.Vec3.IsParallel([vector1.X, vector1.Y, vector1.Z], [
        vector2.X, vector2.Y, vector2.Z],  tolerance)
    if result == 1:
        return True
    if result == 0:
        return False
    assert False, "Invalid return"


def create_orientation_from_z_vector(vector):
    """
    """
    if isinstance(vector, NXOpen.Vector3d):
        return THE_UF_SESSION.Mtx3.InitializeZ([vector.X, vector.Y, vector.Z])
    return THE_UF_SESSION.Mtx3.InitializeZ(vector)


def convert_to_matrix(matrix_array) -> NXOpen.Matrix3x3:
    matrix = NXOpen.Matrix3x3()
    matrix.Xx = matrix_array[0]
    matrix.Xy = matrix_array[1]
    matrix.Xz = matrix_array[2]
    matrix.Yx = matrix_array[3]
    matrix.Yy = matrix_array[4]
    matrix.Yz = matrix_array[5]
    matrix.Zx = matrix_array[6]
    matrix.Zy = matrix_array[7]
    matrix.Zz = matrix_array[8]
    return matrix


def abs_vector(vector) -> NXOpen.Vector3d:
    x = abs(vector.X)
    y = abs(vector.Y)
    z = abs(vector.Z)
    return NXOpen.Vector3d(x, y, z)


def create_note0(detail_number, target_face, tool_face):
    edge_positions = get_edge_positions(tool_face)
    x, y, z = sum_edge_positions(edge_positions)
    length = len(edge_positions)
    origin_in_absolute = NXOpen.Point3d(x / length, y / length, z / length)
    data = target_face.OwningComponent.GetPosition()
    display_part().WCS.SetOriginAndMatrix(data[0], data[1])
    origin = nxopen_point3d_to_list(origin_in_absolute)
    origin_in_target = THE_UF_SESSION.Csys.MapPoint(
        UF_CSYS_ROOT_COORDS, origin, UF_CSYS_ROOT_WCS_COORDS)
    face_vector = normal_vector(target_face)
    temp_orientation = create_orientation_from_z_vector(face_vector)
    orientation = convert_to_matrix(temp_orientation)
    display_part(target_face.Prototype.OwningPart)
    text = create_text_feature(detail_number, NXOpen.Point3d(
        origin_in_target[0], origin_in_target[1], origin_in_target[2]), orientation)
    splines = get_splines(text.GetEntities())
    remove_parameters(work_part(), splines)
    modify_display(splines)
    refset = get_or_create_refset(
        display_part(), REFERENCE_SET)
    refset.AddObjectsToReferenceSet(splines)
    refset.AddObjectsToReferenceSet([target_face.Prototype.GetBody()])


def get_or_create_refset(part, refset_name) -> NXOpen.ReferenceSet:
    for refset in part.GetAllReferenceSets():
        if refset.Name == refset_name:
            return refset
    refset = part.CreateReferenceSet()
    refset.SetName(refset_name)
    return refset


def get_edge_positions(tool_face: NXOpen.Face):
    edge_positions = []
    for edge in tool_face.GetEdges():
        start = start_point_edge(edge)
        end = end_point_edge(edge)
        edge_positions.append(start)
        edge_positions.append(end)
    return edge_positions


def sum_edge_positions(edge_positions):
    x, y, z = 0.0, 0.0, 0.0
    for pos in edge_positions:
        x += pos .X
        y += pos .Y
        z += pos .Z
    return x, y, z


def create_text_feature(detail_number, origin, orientation) -> NXOpen.Features.Text:
    NULL = NXOpen.Features.Text.Null
    text_builder = work_part().Features.CreateTextBuilder(NULL)
    try:
        text_builder.PlanarFrame.AnchorLocation = NXOpen.GeometricUtilities.RectangularFrameBuilder.AnchorLocationType.MiddleCenter
        text_builder.PlanarFrame.WScale = 100.0
        text_builder.PlanarFrame.Length.RightHandSide = str(LENGTH)
        text_builder.PlanarFrame.Height.RightHandSide = str(HEIGHT)
        text_builder.SelectFont(FONT, SCRIPT)
        text_builder.TextString = detail_number
        point2 = work_part().Points.CreatePoint(origin)
        csys = work_part().CoordinateSystems.CreateCoordinateSystem(
            origin, orientation, True)
        text_builder.PlanarFrame.CoordinateSystem = csys
        text_builder.PlanarFrame.UpdateOnCoordinateSystem()
        work_view = work_part().ModelingViews.WorkView
        text_builder.PlanarFrame.AnchorLocator.SetValue(
            point2, work_view, origin)
        text = text_builder.Commit()
    finally:
        text_builder.Destroy()
    return text


def get_splines(entities):
    splines = []
    for t in entities:
        if isinstance(t, NXOpen.Spline):
            splines.append(t)
    return splines


def modify_display(displayable_objects):
    disp_mod = THE_SESSION.DisplayManager.NewDisplayModification()
    try:
        disp_mod.ApplyToAllFaces = True
        disp_mod.ApplyToOwningParts = False
        disp_mod.NewColor = COLOR
        disp_mod.NewLayer = LAYER
        disp_mod.Apply(displayable_objects)
    finally:
        disp_mod.Dispose()


def remove_parameters(part, objects):
    remove_parameters = part.Features.CreateRemoveParametersBuilder()
    try:
        remove_parameters.Objects.Add(objects)
        remove_parameters.Commit()
    finally:
        remove_parameters.Destroy()


def face_is_planar(face: NXOpen.Face) -> bool:
    return face.SolidFaceType == NXOpen.FaceFaceType.Planar


def face_pointing_up_or_down(face: NXOpen.Face) -> bool:
    normal = normal_vector(face)
    unit = unitize_vector(normal)
    abs_vec = abs_vector(unit)
    z_vector = NXOpen.Vector3d(0.0, 0.0, 1.0)
    return is_parallel(abs_vec, z_vector)


def find_interfering_faces(body1: NXOpen.Body, body2: NXOpen.Body) -> list:
    interfering_faces = []
    obj = work_part().AnalysisManager.CreateSimpleInterferenceObject()
    obj.FaceInterferenceType = NXOpen.GeometricAnalysis.SimpleInterference.FaceInterferenceMethod.AllPairs
    obj.InterferenceType = NXOpen.GeometricAnalysis.SimpleInterference.InterferenceMethod.InterferingFaces
    obj.FirstBody.Value = body1
    obj.SecondBody.Value = body2
    obj.PerformCheck()
    result = obj.GetInterferenceResults()
    length = len(result)
    for i in range(0,  length, 2):
        interfering_faces.append((result[i], result[i + 1]))
    obj.Reset()
    return interfering_faces


def touching(face1: NXOpen.Face, face2: NXOpen.Face) -> bool:
    result = THE_UF_SESSION.Modeling.AskMinimumDist(
        face1.Tag, face2.Tag, 0, [0.0, 0.0, 0.0], 0, [0.0, 0.0, 0.0])
    return result[0] == 0.0


def bodies_touch(body1, body2) -> bool:
    result = THE_UF_SESSION.Modeling.CheckInterference(
        body1.Tag, 1, [body2.Tag])
    return result[0] == 3


def main():
    result_target_body = select_solid_body()[1]
    result_tool_bodies = select_solid_bodies()[1]
    t0 = time.time()
    try:
        valid_target_faces = []
        for face in result_target_body.GetFaces():
            if face_is_planar(face) and face_pointing_up_or_down(face):
                valid_target_faces.append(face)
        touching_faces = []
        index = 1
        for tool_body in result_tool_bodies:
            owning_component = tool_body.OwningComponent
            THE_UF_SESSION.Ui.SetPrompt(
                "Processing " + str(index) + " of " + str(len(result_tool_bodies)) + ": " + owning_component.DisplayName)
            index += 1
            if not bodies_touch(result_target_body, tool_body):
                write_line("Unable to find touching faces for body in Component: " +
                           tool_body.OwningComponent.DisplayName)
                continue
            note_created = False
            interfering_faces = find_interfering_faces(
                result_target_body, tool_body)

            target_face = None
            tool_face = None

            for pair in interfering_faces:
                face1 = pair[0]
                face2 = pair[1]
                if not face_is_planar(face1) or not face_is_planar(face2):
                    continue
                if not face_pointing_up_or_down(face1) or not face_pointing_up_or_down(face2):
                    continue
                target_face = pair[0]
                tool_face = pair[1]

            if target_face != None and tool_face != None:
                touching_faces.append((target_face, tool_face))
                note_created = True

            if note_created:
                continue

            for tool_face in tool_body.GetFaces():
                if note_created:
                    break
                if not face_is_planar(tool_face) or not face_pointing_up_or_down(tool_face):
                    continue
                for target_face in valid_target_faces:
                    if note_created:
                        break
                    if touching(target_face, tool_face):
                        touching_faces.append((target_face, tool_face))
                        note_created = True

        for pair in touching_faces:
            target_face = pair[0]
            tool_face = pair[1]
            detail = get_detail_name(tool_face.OwningComponent)
            create_note0(detail, target_face, tool_face)
            tool_face.OwningComponent.Blank()
            write_line("Created note: " + detail)
        return
    finally:
        t1 = time.time()
        write_line("Seconds: " + str(round(t1 - t0, 2)))


if __name__ == '__main__':
    original_display = display_part()
    try:
        main()
    finally:
        display_part(original_display)
