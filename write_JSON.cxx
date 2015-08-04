/* $RCSfile: $
* $Revision: $ $Date: $
* Auth: Nicholas Fay (fayn@rpi.edu)
*
* Copyright (c) 1991-2015 by STEP Tools Inc.
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#include <stp_schema.h>
#include <stix.h>
#include <stixmesh.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <ctype.h>

#include "stp2webgl.h"
#include "json.h"

using namespace Json;
using namespace std;

// transfor moved into stix in latest version
#ifndef LATEST_STDEV
#define stix_get_transform stixmesh_get_transform
#endif

// This function writes a lightweight JSON description of the STEP
// product structure and faceted shapes using the format described
// below.
//
// http://www.steptools.com/support/stdev_docs/stixmesh/XMLFormat.html
//
// The product structure is written in a first pass and all of the
// shapes are submitted to the facetter.  Then the shell facets are
// written as they become available and then the facetted data is
// released.  This is a more complex arrangement than just facetting
// everything in one batch, but it may be more memory efficient.
//


extern int write_JSON(stp2webgl_opts * opts);

//======================================================================

static string append_ref(RoseObject * obj)
{
	if (obj->isa(ROSE_DOMAIN(RoseUnion))) {
		obj = rose_get_nested_object(ROSE_CAST(RoseUnion, obj));
	}

	if (!obj->entity_id()) {
		printf("No entity id for %p: %s\n", obj, obj->domain()->name());
		exit(2);
	}

	char numstr[21];
	string id = "id";
	string full_id = id + itoa(obj->entity_id(), numstr, 10);


	return full_id;
}


static FILE * open_dir_file(const char * dir, const char * fname)
{
	RoseStringObject path = dir;
	path.cat("/");
	path.cat(fname);

	return fopen(path, "w");
}


//======================================================================
// STEP PMI Annotations -- GD&T and construction planes
// 



static void write_poly_point(double vals[3], Value &poly_point)
{
	Value polyline;
	char buff[64];
	sprintf(buff, "%.15g %.15g %.15g", vals[0], vals[1], vals[2]);
	polyline["l"] = buff;
	poly_point["p"] = polyline;
}

static void write_poly_point(ListOfDouble * vals, Value &poly_point)
{
	Value polyline;
	char buff[64];
	sprintf(buff, "%.15g %.15g %.15g", vals->get(0), vals->get(1), vals->get(2));
	polyline["l"] = buff;
	poly_point["p"] = polyline;
}


static void append_step_curve(stp_curve * c, Value &stepArr)
{
	/* May also want to handle composite curves here in this routine, since
	* a composite curve may contain via points, and thus run into trouble.*/
	unsigned i, sz;

	if (c->isa(ROSE_DOMAIN(stp_polyline))) {

		stp_polyline * poly = ROSE_CAST(stp_polyline, c);
		ListOfstp_cartesian_point * pts = poly->points();

		if (!pts || pts->size() < 2)
			return;

		Value plines;
		for (i = 0, sz = pts->size(); i<sz; i++) {
			stp_cartesian_point * pt = pts->get(i);
			ListOfDouble * vals = pt->coordinates();
			Value pline;
			write_poly_point(vals,pline);
			plines.append(pline);
		}
		stepArr["polyline"] = plines;
	}

}

static void append_annotation(
	stp_geometric_set * gset,
	Value &annot
	)
{
	SetOfstp_geometric_set_select * elems = gset->elements();
	if (!elems) return;

	unsigned i, sz;
	Value stepArr;
	for (i = 0, sz = elems->size(); i<sz; i++) {
		stp_geometric_set_select * sel = elems->get(i);
		Value annot2;
		if (sel->is_curve()) {
			append_step_curve(sel->_curve(), annot2);
			stepArr.append(annot2);
		}
		else if (sel->is_point()) {
		}
		else if (sel->is_surface()) {
		}

	}
	annot.append(stepArr);

}

static void append_annotation(
	stp_representation_item * it,
	Value &annot
	)
{
	if (!it) return;

	if (it->isa(ROSE_DOMAIN(stp_annotation_plane)))
	{
		// do we do anything special with the plane?
		unsigned i, sz;
		stp_annotation_plane * ap = ROSE_CAST(stp_annotation_plane, it);

		Value annotArr;
		for (i = 0, sz = ap->elements()->size(); i<sz; i++)
		{
			// either draughting_callout or styled_item, both of which
			// are rep items.
			stp_representation_item * elem =
				ROSE_CAST(stp_representation_item,
				rose_get_nested_object(ap->elements()->get(i)));
			Value annot2;
			append_annotation(elem,annot2);
			annotArr.append(annot2);
		}
		annot.append(annotArr);
	}

	else if (it->isa(ROSE_DOMAIN(stp_annotation_occurrence)))
	{
		stp_annotation_occurrence * ao
			= ROSE_CAST(stp_annotation_occurrence, it);

		if (!ao->item()) return;

		if (ao->item()->isa(ROSE_DOMAIN(stp_geometric_set))) {
			Value annotation2;
			append_annotation(ROSE_CAST(stp_geometric_set, ao->item()), annotation2);
			annot.append(annotation2);
		}
	}

	else {
		printf("append_annotation unimplemented case: %s\n",
			it->domain()->name());
	}

}


static void append_model_body(
	stp2webgl_opts * opts,
	stp_representation * model,
	Value &model_bod
	)
{	
	string id = append_ref(model);
	model_bod["id"] = id;

	SetOfstp_representation_item * items = model->items();
	unsigned sz = items->size();
	Value rep_items;
	for (unsigned i = 0; i<sz; i++) {
		Value rep_item;
		stp_representation_item * it = items->get(i);
		append_annotation(it,rep_item);
		rep_items.append(rep_item);
	}
	model_bod["annotation"] = rep_items;
}


static void append_model(
	stp2webgl_opts * opts,
	stp_representation * model,
	Value &model2
	)
{
	if (!model || rose_is_marked(model)) return;
	rose_mark_set(model);

	if (!opts->do_split)
		append_model_body(opts, model, model2);
	else
	{
		Value tempVal;
		char fname[100];
		sprintf(fname, "annotation_id%lu.JSON", model->entity_id());

		string id = append_ref(model);

		FILE * fd = open_dir_file(opts->dstdir, fname);

		append_model_body(opts, model, tempVal);
		model2 = tempVal;
		tempVal["id"] = id;
		tempVal["href"] = fname;

		Value exportVal;
		exportVal["annotation"] = tempVal;

		string output = exportVal.toStyledString();
		const char * c = output.c_str();
		fputs(c, fd);

		fclose(fd);
	}
}


static void append_step_curve(
	stp_representation * rep,
	stp_bounded_curve * curve,
	Value &curves
	)
{
	StixMeshNurbs nurbs;
	stixmesh_create_bounded_curve(&nurbs, curve, rep);

	RoseBoundingBox bbox;

	if (!nurbs.getConvexHull(&bbox)) {
		printf("Could not get convec hull of curve, skipping");
		return;
	}
	double tol = bbox.diagonal() / 100.;

	rose_real_vector u_vals;
	nurbs.extractTolerancedPoints(&u_vals, tol, 1);

	unsigned i, sz;
	Value polyArr;
	for (i = 0, sz = u_vals.size(); i<sz; i++) {
		double xyz[3];
		Value poly_point;
		nurbs.eval(xyz, u_vals[i]);
		write_poly_point(xyz, poly_point);
		polyArr.append(poly_point);
	}
	curves["polyline"] = polyArr;
}

static void append_rep_item(
	stp_representation * rep,
	stp_representation_item * it,
	Value &item
	)
{
	if (!it)
		return;

	if (it->isa(ROSE_DOMAIN(stp_bounded_curve))) {
		append_step_curve(rep, ROSE_CAST(stp_bounded_curve, it), item);
	}
}

static void append_constructive_geom_body(
	stp2webgl_opts * opts,
	stp_constructive_geometry_representation * cgr,
	Value &geom
	)
{
	unsigned i, sz;
	string id = append_ref(cgr);

	Value repArr;

	SetOfstp_representation_item * items = cgr->items();
	for (i = 0, sz = items->size(); i<sz; i++) {
		stp_representation_item * it = items->get(i);
		Value item;
		append_rep_item(cgr, it, item);
		repArr.append(item);
	}
	geom["annotation"] = repArr;
}


static void append_constructive_geom(
	stp2webgl_opts * opts,
	stp_constructive_geometry_representation * cg,
	Value &con_geom
	)
{
	Value cgeometry;
	if (!cg || rose_is_marked(cg)) return;
	rose_mark_set(cg);

	if (!opts->do_split){
		append_constructive_geom_body(opts, cg, cgeometry);
		con_geom.append(cgeometry);
	}
	else
	{
		Value exportVal;
		Value con_geom_body;

		char fname[100];
		sprintf(fname, "constructive_id%lu.JSON", cg->entity_id());
		string id = append_ref(cg);

		FILE * fd = fopen(fname, "w");

		append_constructive_geom_body(opts, cg, con_geom_body);

		exportVal["annotation"] = con_geom_body;
		con_geom["annotation"] = con_geom_body;

		con_geom_body["id"] = id;
		con_geom_body["href"] = fname;

		string output = exportVal.toStyledString();
		const char * c = output.c_str();
		fputs(c, fd);

		fclose(fd);
	}
}



static void append_annotation_refs(
	stp_representation * rep,
	Value &annot
	)
{
	/* FIXME - this generates both constructive geometry and
	* annotations.  theses should be split out, but are are not
	* (currently) doing so since that would require updating the
	* webgl javascript.
	*/
	unsigned i, sz;
	unsigned cnt = 0;

	if (!rep->isa(ROSE_DOMAIN(stp_shape_representation)))
		return;

	StixMeshRepresentationVec * models = stixmesh_get_draughting_models(
		ROSE_CAST(stp_shape_representation, rep)
		);

	StixMeshConstructiveGeomVec * cgeom =
		stixmesh_get_constructive_geometry(rep);

	if (!models && !cgeom)
		return;

	string id_chain;
	if (models) {
		for (i = 0, sz = models->size(); i<sz; i++) {
			if (cnt){
				string x = append_ref(models->get(i));
				cnt++;
				id_chain = id_chain + " " + x;
			}
		}
	}

	if (cgeom) {
		for (i = 0, sz = cgeom->size(); i<sz; i++) {
			if (cnt){
				string y = append_ref(cgeom->get(i));
				cnt++;
				id_chain = id_chain + " " + y;
			}
		}
	}
	annot["annotation"] = id_chain;
}



static void append_annotations(
	stp2webgl_opts * opts,
	stp_representation * rep,
	Value &models2,
	Value &cgeom2
	)
{
	unsigned i, sz;

	if (!rep->isa(ROSE_DOMAIN(stp_shape_representation)))
		return;

	StixMeshRepresentationVec * models = stixmesh_get_draughting_models(
		ROSE_CAST(stp_shape_representation, rep));

	StixMeshConstructiveGeomVec * cgeom =
		stixmesh_get_constructive_geometry(rep);

	if (!models && !cgeom)
		return;

	if (models) {
		Value modelArr;
		for (i = 0, sz = models->size(); i<sz; i++) {
			stp_representation * model = models->get(i);
			Value model2;
			append_model(opts, model, model2);
			modelArr.append(model2);
		}
		models2["m"] = modelArr;
	}

	if (cgeom) {
		Value cgeomArr;
		for (i = 0, sz = cgeom->size(); i<sz; i++) {
			stp_constructive_geometry_representation * cg = cgeom->get(i);
			Value cgeomz;
			append_constructive_geom(opts, cg, cgeomz);
			cgeomArr.append(cgeomz);
		}
		cgeom2["cg"] = cgeomArr;
	}

}


//======================================================================
// Queue STEP Representation -- Add shape data to facetter queue and
// write forward information about the step shell.
//
void append_shell_refs(
	stp_representation * rep,
	Value &parent
	)
{
	unsigned i, sz;
	unsigned count = 0;
	SetOfstp_representation_item * items = rep->items();
	string shell_ids;

	for (i = 0, sz = items->size(); i<sz; i++)
	{
		stp_representation_item * it = items->get(i);
		if (!StixMeshStpBuilder::isShell(rep, it)){
			count++;
			continue;
		}

		string shell_id = append_ref(it);
		if ((i-count) == 0){
		    shell_ids = shell_id;
		}
		else{
		    shell_ids = shell_ids + " " + shell_id;
		}
	}
	if ((items->size() - count) > 0){
		parent["shell"] = shell_ids;
	}
}


static void append_asm_child(
	stp2webgl_opts * opts,
	RoseObject * rel,
	Value &parent
	)
{
	StixMgrAsmRelation * mgr = StixMgrAsmRelation::find(rel);
	if (!mgr) return;

	stp_representation * child = mgr->child;

	string ref = append_ref(child);
	Value temp;
	temp["ref"] = ref;

	unsigned i;
	RoseXform xform = stix_get_transform(mgr);
	string xformString;
	for (i = 0; i<16; i++) {
	    if (i == 16){
		char buff[64];
		sprintf(buff, "%.15g", xform.m[i]);
		xformString = xformString + buff;
	    }
	    else{
		char buff[64];
		sprintf(buff, "%.15g ", xform.m[i]);
		xformString = xformString + buff;
	    }
	}

	temp["xform"] = xformString;
	parent["child"].append(temp);
}


void queue_shapes(
	stp2webgl_opts * opts,
	StixMeshStpAsyncMaker * mesher,
	stp_representation * rep,
	Value &shapeArray
	)
{
	unsigned i, j, sz;

	if (!rep || rose_is_marked(rep)) return;
	rose_mark_set(rep);

	StixMgrAsmShapeRep * mgr = StixMgrAsmShapeRep::find(rep);
	if (!mgr) return;

	// Write facets by default, unless we have a list of reps.  In the
	// latter case only write facets if a rep is in the list.

	int do_facets = 1;
	if (opts->root_ids.size())
	{
		unsigned eid = rep->entity_id();
		do_facets = 0;

		for (i = 0, sz = opts->root_ids.size(); i<sz; i++)
		{
			if (opts->root_ids[i] == eid) {
				do_facets = 1;
				break;
			}
		}
	}

	Value shapeVal;

	string id = append_ref(rep);
	shapeVal["id"] = id;

	RoseUnit unit = stix_get_context_length_unit(rep);
	if (unit != roseunit_unknown)
	{
		char buff[20];
		sprintf(buff, "%s %f",stix_get_unit_name(unit), rose_get_measure_as_unit(1., unit, roseunit_m));
		shapeVal["unit"] = buff;
	}

	if (do_facets){
		append_shell_refs(rep, shapeVal);
	}

	append_annotation_refs(rep,shapeVal);

	for (j = 0, sz = mgr->child_rels.size(); j<sz; j++)
		append_asm_child(opts, mgr->child_rels[j], shapeVal);

	for (j = 0, sz = mgr->child_mapped_items.size(); j<sz; j++)
		append_asm_child(opts, mgr->child_mapped_items[j], shapeVal);

	shapeArray.append(shapeVal);

	for (j = 0, sz = mgr->child_rels.size(); j<sz; j++) {
		StixMgrAsmRelation * rm = StixMgrAsmRelation::find(
			mgr->child_rels[j]
			);
		if (rm) queue_shapes(opts, mesher, rm->child, shapeArray);
	}

	for (j = 0, sz = mgr->child_mapped_items.size(); j<sz; j++) {
		StixMgrAsmRelation * rm = StixMgrAsmRelation::find(
			mgr->child_mapped_items[j]
			);
		if (rm) queue_shapes(opts, mesher, rm->child, shapeArray);
	}


	if (do_facets)
	{
		SetOfstp_representation_item * items = rep->items();
		unsigned sz = items->size();

		for (unsigned i = 0; i<sz; i++) {
			stp_representation_item * ri = items->get(i);

			if (StixMeshStpBuilder::canMake(rep, ri))
				mesher->startMesh(rep, ri, &opts->mesh);
		}
		Value shap;
		Value shap2;
		append_annotations(opts, rep, shap,shap2);
		if (!shap.isNull()){
			shapeArray.append(shap);
		}
		if (!shap2.isNull()){
			shapeArray.append(shap2);
		}
	}
}



//======================================================================
// Write Triangle Data for a Facetted Shell
//
static void append_facet(
	const RoseMeshFacetSet * fs,
	unsigned fidx,
	int write_normal,
	Value &facets
	)
{
	const RoseMeshFacet * f = fs->getFacet(fidx);
	if (!f) return;

	char buff[64];
	sprintf(buff, "%d, %d, %d", f->verts[0], f->verts[1], f->verts[2]);
	facets["v"] = buff;

	if (write_normal) {
		// facet_normal_now_computed_in_latest_versions
#ifdef LATEST_STDEV
		double fnorm[3];
		fs->getFacetNormal(fnorm, fidx);
#else
		const double * fnorm = fs->getNormal(f->facet_normal);
#endif
		char fnorms[64];
		sprintf(fnorms, "%.15g, %.15g, %.15g", fnorm[0], fnorm[1], fnorm[2]);
		facets["fn"] = fnorms;
	}

	Value normals;
	for (unsigned j = 0; j<3; j++) {
		// facet_normal_now_computed_in_latest_versions
#ifdef LATEST_STDEV
		const double * normal = fs->getNormal(f->normals[j]);
#else
		const double * normal = fs->getNormal(f->vert_normals[j]);
#endif
		if (normal)
		{
			Value sing_norm;
			char buff2[256];
			sprintf(buff2, "%.15g %.15g %.15g", normal[0], normal[1], normal[2]);
			sing_norm["d"] = buff2;
			normals.append(sing_norm);
		}
	}
		facets["n"] = normals;
}


void append_shell_facets(
	const StixMeshStp * shell,
	Value &shellVal
	)
{
	int WRITE_NORMAL = 0;
	unsigned i, sz;
	unsigned j, szz;
	const RoseMeshFacetSet * facets = shell->getFacetSet();


	string id = append_ref(shell->getStepSolid());
	shellVal["id"] = id;

	unsigned dflt_color = stixmesh_get_color(shell->getStepSolid());
	if (dflt_color != ROSE_MESH_NULL_COLOR){
		char buff[] = "rrggbb ";
		sprintf(buff, "%06x", dflt_color);
		shellVal["color"] = buff;
	}

	Value verts;
	for (i = 0, sz = facets->getVertexCount(); i<sz; i++)
	{
		const double * pt = facets->getVertex(i);
		Value point;
		char buff[64];
		sprintf(buff, "%.15g %.15g %.15g", pt[0], pt[1], pt[2]);
		point["p"] = buff;
		verts.append(point);
	}

	Value vert_arr;
	vert_arr["v"] = verts;
	shellVal["verts"] = vert_arr;

	// The facet set has all of the facets for the shell.  Break it up
	// into groups by step face.

	Value facets2;
	string color2;
	for (i = 0, sz = shell->getFaceCount(); i<sz; i++)
	{
		const StixMeshStpFace * fi = shell->getFaceInfo(i);
		unsigned first = fi->getFirstFacet();
		unsigned color = stixmesh_get_color(fi->getFace());
		if (first == ROSE_NOTFOUND)
			continue;

		// Always tag the face with a color, unless everything is null. 
		if (color == ROSE_MESH_NULL_COLOR)
			color = dflt_color;

		if (color != ROSE_MESH_NULL_COLOR){
			char buff2[] = "rrggbb ";
			sprintf(buff2, "%06x", dflt_color);
			color2 = buff2;
			facets2["color"] = color2;
		}
		Value tempfacet;
		for (j = 0, szz = fi->getFacetCount(); j<szz; j++) {
			append_facet(facets, j + first, WRITE_NORMAL, tempfacet);
			facets2["f"].append(tempfacet);
		}
		shellVal["facets"].append(facets2);
	}
}


static void export_shell(
	stp2webgl_opts * opts,
	const StixMeshStp * shell,
	Value &shell_arr
	)
{
	if (!shell) return;

	Value shellVal;
	Value exportVal;
	Value exportVal2;

	if (!opts->do_split) {
		append_shell_facets(shell,shellVal);
		shell_arr.append(shellVal);
	}
	else
	{
		unsigned i, sz;
		RoseBoundingBox bbox;
		const RoseMeshFacetSet * facets = shell->getFacetSet();

		// compute the bounding box for the shell
		for (i = 0, sz = facets->getVertexCount(); i<sz; i++)
		{
			const double * pt = facets->getVertex(i);
			bbox.update(pt);
		}

		char fname[100];
		sprintf(fname, "shell_id%lu.JSON", shell->getStepSolid()->entity_id());

		shellVal["id"] = append_ref(shell->getStepSolid());
		shellVal["size"] = facets->getFacetCount();

		Value bboxArr;
		//Turn the 6 bbox values into an array and set equal to bbox
		bboxArr.append(bbox.minx());
		bboxArr.append(bbox.miny());
		bboxArr.append(bbox.minz());
		bboxArr.append(bbox.maxx());
		bboxArr.append(bbox.maxy());
		bboxArr.append(bbox.maxz());

		shellVal["bbox"] = bboxArr;

		// append the area 
		double area = 0.;
		for (i = 0, sz = shell->getFaceCount(); i<sz; i++) {
			const StixMeshStpFace * face = shell->getFaceInfo(i);
			area += face->getArea();
		}

		shellVal["area"] = area;
		shellVal["href"] = fname;

		shell_arr.append(shellVal);

		FILE * fd = open_dir_file(opts->dstdir, fname);

		append_shell_facets(shell,exportVal);
		exportVal2["shell"] = exportVal;

		FastWriter writer;

		string output = writer.write(exportVal2);
		const char * c = output.c_str();
		fputs(c, fd);

		fclose(fd);
	}
}




//======================================================================
// Write STEP Product Structure
//

static const char * get_product_name(stp_next_assembly_usage_occurrence* nauo)
{
	stp_product_definition * pd = stix_get_related_pdef(nauo);
	stp_product_definition_formation * pdf = pd ? pd->formation() : 0;
	stp_product * prod = pdf ? pdf->of_product() : 0;

	return prod ? prod->name() : 0;
}

int nauo_product_cm(const void* a, const void* b)
{
	const char * name_a = get_product_name(
		*(stp_next_assembly_usage_occurrence**)a
		);

	const char * name_b = get_product_name(
		*(stp_next_assembly_usage_occurrence**)b
		);

	if (name_a == name_b) return 0;
	if (!name_a) return +1;
	if (!name_b) return -1;

	return strcmp(name_a, name_b);
}

static void append_stplink(
	stp2webgl_opts * opts,
	stp_product_definition * pd,
	Value &obj
	)
{
	// write reference to a component stepfile.  Only used when
	// splitting a large assembly into a a collection of part files
	// and when a file exists for this particular product.
	// 
	if (!opts->do_split) return;

	RoseStringObject name("part");

	stp_product_definition_formation * pdf = pd->formation();
	stp_product * p = pdf ? pdf->of_product() : 0;

	char * pname = p ? p->name() : 0;
	if (!pname || !*pname) pname = (char *) "none";

	if (!name.is_empty()) name += "_";
	name += pname;

	// change whitespace and other non filesystem safe
	// characters to underscores
	//
	char * c = name;
	while (*c) {
		if (isspace(*c)) *c = '_';
		if (*c == '?')   *c = '_';
		if (*c == '/')   *c = '_';
		if (*c == '\\')  *c = '_';
		if (*c == ':')   *c = '_';
		if (*c == '"')   *c = '_';
		if (*c == '\'')  *c = '_';
		c++;
	}

	if (pd->design()->fileExtension()) {
		name += ".";
		name += pd->design()->fileExtension();
	}

	RoseStringObject path = opts->dstdir;
	path += "/";
	path += name;

	// expects the component files to be present already.  Do not
	// generate link if not there.
	if (rose_file_exists(path)) {
		obj["step"] = name.as_char();
	}
	else{
	    obj["step"] = name.as_char();
	}
}


static void export_product(
	stp2webgl_opts * opts,
	stp_product_definition * pd,
	Value &temparray
	)
{
	unsigned i, sz;

	if (!pd || rose_is_marked(pd)) return;
	rose_mark_set(pd);

	StixMgrAsmProduct * mgr = StixMgrAsmProduct::find(pd);

	string id = append_ref(pd);
	Value tempObj;
	tempObj["id"] = id;

	append_stplink(opts, pd, tempObj);

	stp_product_definition_formation * pdf = pd->formation();
	stp_product * prod = pdf->of_product();

	tempObj["name"] = prod->name();


	if (mgr->shapes.size()) {
		string temp;
		for (i = 0, sz = mgr->shapes.size(); i<sz; i++)
		{
			string shp_id = append_ref(mgr->shapes[i]);
			if (i == 0){
			    temp = shp_id;
			}
			else{
			    temp = temp + " " + shp_id;
			}
		}
		tempObj["shape"] = temp;
	}

	if (mgr->child_nauos.size())
	{
		qsort(mgr->child_nauos._buffer(),
			mgr->child_nauos.size(),
			sizeof(stp_next_assembly_usage_occurrence*),
			&nauo_product_cm
			);
		string id_string;
		for (i = 0, sz = mgr->child_nauos.size(); i<sz; i++)
		{
			stp_next_assembly_usage_occurrence * nauo = mgr->child_nauos[i];
			string id2 = append_ref(stix_get_related_pdef(mgr->child_nauos[i]));
			if (i == 0){
			    id_string = id2;
			}
			else{
			    id_string = id_string + " " + id2;
			}
		}
		tempObj["children"] = id_string;
	}

	temparray.append(tempObj);

	for (i = 0, sz = mgr->child_nauos.size(); i<sz; i++)
	{
		stp_next_assembly_usage_occurrence * nauo = mgr->child_nauos[i];
		export_product(opts, stix_get_related_pdef(nauo), temparray);
	}
}


// ======================================================================


int write_JSON(stp2webgl_opts * opts)
{
	FILE * JSONout = 0;
	RoseStringObject index_file;
	unsigned i, sz;

	if (opts->do_split)
	{
		opts->dstdir = opts->dstfile;
		if (!opts->dstdir)
			opts->dstdir = "step_data";

		index_file = opts->dstdir;
		index_file.cat("/index.JSON");
		opts->dstfile = index_file;

		if (!rose_dir_exists(opts->dstdir) &&
			(rose_mkdir(opts->dstdir) != 0)) {
			printf("Cannot create directory %s\n", opts->dstdir);
			return 2;
		}
	}

	if (opts->dstfile)
	{
		JSONout = rose_fopen(opts->dstfile, "w");
		if (!JSONout) {
			printf("Could not open output file\n");
			return 2;
		}
	}

	if (JSONout == 0) {
		JSONout = fopen("index.JSON", "w");
	}

	// The RoseOutputFile class is a data stream
	// class that the JSON file writes to. 
	//
	RoseOutputFile JSONfile(JSONout, opts->dstfile ? opts->dstfile : "JSON file");
	Value rootval;
	Value rootval2;

	for (i = 0, sz = opts->root_prods.size(); i<sz; i++)
	{
		stp_product_definition * pd = opts->root_prods[i];	
		string val2 = append_ref(pd);
		rootval2["-root"] = val2;
	}


	rose_mark_begin();
	
	Value temparray = arrayValue;
	for (i = 0, sz = opts->root_prods.size(); i<sz; i++)
	{
		export_product(opts, opts->root_prods[i], temparray);
	}
	rootval2["product"] = temparray;

	// Schedule each solid for faceting, which will happen in child
	// threads and then write each shell as it becomes available.
	StixMeshStpAsyncMaker mesher;
	StixMeshStp * mesh;

	Value temparray2 = arrayValue;
	for (i = 0, sz = opts->root_prods.size(); i<sz; i++)
	{
		unsigned j, szz;
		StixMgrAsmProduct * mgr = StixMgrAsmProduct::find(
			opts->root_prods[i]
			);

		for (j = 0, szz = mgr->shapes.size(); j<szz; j++) {
			queue_shapes(opts, &mesher, mgr->shapes[j], temparray2);
		}
	}
	rootval2["shape"] = temparray2;

	Value temparray3 = arrayValue;
	
	while ((mesh = mesher.getResult(1)) != 0)
	{
		export_shell(opts, mesh, temparray3);
		delete mesh;
	}
	rootval2["shell"] = temparray3;

	rootval["step_assembly"] = rootval2;

	//string output = rootval.toStyledString();
	FastWriter writer;

	string output = writer.write(rootval);
	const char * c = output.c_str();
	fputs(c,JSONout);

	fclose(JSONout);

	JSONfile.flush();
	rose_mark_end();
	return 0;
}