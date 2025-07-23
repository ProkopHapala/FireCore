#version 330 core
out vec4 FragColor;

in vec4 sphere_obj_world; // Sphere center (xyz) and radius (w) in world space
in vec3 fpos_world;       // Fragment position on the bounding mesh in world space (ray endpoint)
in vec4 fColor;           // Color from vertex shader

// Uniforms
uniform vec3 viewPos;      // Camera position in world space (ray origin)
uniform mat4 projection;   // For gl_FragDepth
uniform mat4 view;         // For gl_FragDepth

vec2 rayPointDist( vec3 ray_origin, vec3 ray_direction_normalized, vec3 point ){
	vec3 to_point  = point - ray_origin;
	float t_closest_approach = dot( ray_direction_normalized, to_point ); // Projection of to_point onto ray
    // If t_closest_approach < 0, the closest point is behind the ray origin,
    // but for rays starting at viewPos and going "forward", this is usually fine.
    vec3 point_on_ray_closest = ray_origin + t_closest_approach * ray_direction_normalized;
	float dist_sq_to_center = dot( point - point_on_ray_closest, point - point_on_ray_closest );
	return vec2( t_closest_approach, sqrt(dist_sq_to_center) ); // t, distance_to_center
}

void main()
{
    vec3 ray_origin    = viewPos;
    vec3 ray_direction = normalize(fpos_world - viewPos);
    vec3  sphere_center_w = sphere_obj_world.xyz;
    float sphere_radius_w = sphere_obj_world.w;

    // res.x = t_closest_approach (distance along ray to point of closest approach)
    // res.y = distance_from_ray_to_sphere_center (shortest distance)
    vec2 res = rayPointDist( ray_origin, ray_direction, sphere_center_w );
    float t_closest = res.x;
    float dist_ray_to_center = res.y;

    if (dist_ray_to_center >= sphere_radius_w) {
        // Ray misses the sphere entirely (closest approach is outside radius)
        discard;
    } else {
        float r_w = dist_ray_to_center/sphere_radius_w;
        //float density = exp(-4.*r_w*r_w);
        float density = exp(-2.*r_w*r_w);

        float alpha = density * fColor.a;

        FragColor = vec4(fColor.rgb, alpha);

        // Depth: Use the depth of the point on the ray closest to the sphere center,
        // but only if it's within the sphere.
        vec3 P_closest_on_ray_world = ray_origin + t_closest * ray_direction;
        vec4 clip_space_pos = projection * view * vec4(P_closest_on_ray_world, 1.0);
        gl_FragDepth = (clip_space_pos.z / clip_space_pos.w * 0.5) + 0.5;
    }
}