#define SOLUTION_CYLINDER_AND_PLANE
#define SOLUTION_MATERIAL
#define SOLUTION_SHADOW
#define SOLUTION_REFLECTION_REFRACTION
#define SOLUTION_FRESNEL

precision highp float;
uniform float time;

struct	PointLight 
{
    vec3	position;
    vec3	color;
};

struct	Material 
{
    vec3	diffuse;
    vec3	specular;
    float	glossiness;
#ifdef SOLUTION_MATERIAL
	float	reflection;
	float	refraction;
	float	IOR;
#endif
};

struct	Sphere 
{
    vec3		position;
    float		radius;
    Material	material;
};

struct	Plane 
{
    vec3		normal;
    float		d;
    Material	material;
};

struct	Cylinder 
{
    vec3		position;
    vec3		direction;  
    float		radius;
    Material	material;
};

const int	lightCount = 2;
const int	sphereCount = 3;
const int	planeCount = 1;
const int	cylinderCount = 2;

struct	Scene 
{
    vec3						ambient;
    PointLight[lightCount]		lights;
    Sphere[sphereCount]			spheres;
    Plane[planeCount]			planes;
    Cylinder[cylinderCount]		cylinders;
};

struct	Ray 
{
    vec3	origin;
    vec3	direction;
};

// Contains all information pertaining to a ray/object intersection
struct	HitInfo 
{
    bool		hit;
    float		t;
    vec3		position;
    vec3		normal;
    Material	material;
};

HitInfo	getEmptyHit() 
{
	return	HitInfo
	(
      	false, 
      	0.0, 
      	vec3(0.0), 
      	vec3(0.0), 
#ifdef SOLUTION_MATERIAL
		// Update the constructor call
		Material(vec3(0.0), vec3(0.0), 0.0, 0.0, 0.0, 0.0)
#else
		Material(vec3(0.0), vec3(0.0), 0.0)
#endif
	);
}

// Sorts the two t values such that t1 is smaller than t2
void sortT(inout float t1, inout float t2) 
{
  	// Make t1 the smaller t
    if (t2 < t1) 
	{
		float temp = t1;
		t1 = t2;
		t2 = temp;
    }
}

// Tests if t is in an interval
bool isTInInterval(const float t, const float tMin, const float tMax) 
{
	return t > tMin && t < tMax;
}

// Get the smallest t in an interval
bool getSmallestTInInterval(float t0, float t1, const float tMin, const float tMax, inout float smallestTInInterval) 
{
	sortT(t0, t1);
	// As t0 is smaller, test this first
	if (isTInInterval(t0, tMin, tMax)) 
	{
		smallestTInInterval = t0;
        return true;
	}
  
	// If t0 was not in the interval, still t1 could be
	if (isTInInterval(t1, tMin, tMax)) 
	{
		smallestTInInterval = t1;
		return true;
	}  
	// None was
	return false;
}

HitInfo intersectSphere(const Ray ray, const Sphere sphere, const float tMin, const float tMax) 
{       
	vec3 to_sphere = ray.origin - sphere.position;
    float a = dot(ray.direction, ray.direction);
    float b = 2.0 * dot(ray.direction, to_sphere);
    float c = dot(to_sphere, to_sphere) - sphere.radius * sphere.radius;
    float D = b * b - 4.0 * a * c;
    if (D > 0.0)
    {
		float t0 = (-b - sqrt(D)) / (2.0 * a);
		float t1 = (-b + sqrt(D)) / (2.0 * a);
      
      	float smallestTInInterval;

      	if (!getSmallestTInInterval(t0, t1, tMin, tMax, smallestTInInterval)) 
		{
			return getEmptyHit();
		}
      
      	vec3 hitPosition = ray.origin + smallestTInInterval * ray.direction;      

      	vec3 normal = 
		  length(ray.origin - sphere.position) < sphere.radius + 0.001 ? 
		  -normalize(hitPosition - sphere.position) : 
		  normalize(hitPosition - sphere.position);      

        return HitInfo
		(
          	true,
          	smallestTInInterval,
          	hitPosition,
          	normal,
          	sphere.material
        );
    }
    return getEmptyHit();

}

HitInfo intersectPlane(const Ray ray,const Plane plane, const float tMin, const float tMax) 
{
#ifdef SOLUTION_CYLINDER_AND_PLANE 
	// ------ Add your plane intersection code here

	// Ray: P = P0 (Original_position) + t * V (Ray_Direction)
	// Plane: p * N (Plane_normal) + d (Plane_position, distance to the coordinate origin)
	// Substituting for P,  (P0 + tV) * N + d = 0
	// So,  t = -(P0 * N + d) / (V * N)
	//		P = P0 + t * V

	// For D less than 0, the plane will occur below the ray origin.
	float D = dot(ray.direction,plane.normal);
	if (D < 0.0)	
	{
		float t = -(dot(ray.origin,plane.normal) + plane.d) / dot(ray.direction,plane.normal);
		vec3 hitPosition = ray.origin + t * ray.direction;
		return HitInfo
		(
		true,
		t,
		hitPosition,
		plane.normal,
		plane.material
		);
	}

    return getEmptyHit();
#endif 

    return getEmptyHit();
}

float lengthSquared(vec3 x) 
{
	return dot(x, x);
}

HitInfo intersectCylinder(const Ray ray, const Cylinder cylinder, const float tMin, const float tMax) 
{
#ifdef SOLUTION_CYLINDER_AND_PLANE 

    // ------ Add your cylinder intersection code here

	// In general, an infinite cylinder along Z-axis has equation as x^2 + y^2 = r^2. For more general scenario, an infinite cylinder along line (Pa + Va * t) can be regarded as :
	// (q - Pa - (Va, q - Pa) * Va)^2 - r^2 = 0. 
	// Where, Pa is the position of cylinder in this coursework, Va is the direction of the line/Cylinder. q is for the position on the cylinder. 
	// And then q can be subsyituted as P + V *t, where P is the position of the ray, and v is the direction of the ray.
	// Then we get, (P + V * t - Pa - (Va, P + V * t - Pa) * Va)^2 - r^2 = 0
	// In this equation, term (P + V * t - Pa) is a vector of the cylinder position to a point on the cylinder; 
	// and term ((Va, P + V * t - Pa) * Va) is a vector of the cylinder position to the point on the axis of the cylinder with a normal direction to the point on the cylinder.
	// Reduce the equation to the form of A * t^2 + B * t + C = 0,
	// a = (V - (V, Va) * Va)^2
	// b = 2 * (V - (V, Va) * Va, P - Pa - (P - Pa, Va) * Va)
	// c = (P - Pa - (P - Pa, Va) * Va)^2 - r^2

    vec3 to_cylinder = ray.origin - cylinder.position;
    float a = dot((ray.direction - dot(ray.direction, cylinder.direction) * cylinder.direction),(ray.direction - dot(ray.direction, cylinder.direction) * cylinder.direction));
    float b = 2.0 * dot((ray.direction - dot(ray.direction, cylinder.direction) * cylinder.direction),to_cylinder - dot(to_cylinder, cylinder.direction) * cylinder.direction);
    float c = dot((to_cylinder - dot(to_cylinder, cylinder.direction) * cylinder.direction),(to_cylinder - dot(to_cylinder, cylinder.direction) * cylinder.direction)) - cylinder.radius * cylinder.radius;
    float D = b * b - 4.0 * a * c; 
	if (D > 0.0)
    {
		float t0 = (-b - sqrt(D)) / (2.0 * a);
		float t1 = (-b + sqrt(D)) / (2.0 * a);
      	float smallestTInInterval;

      	if (!getSmallestTInInterval(t0, t1, tMin, tMax, smallestTInInterval)) 
		{
			return getEmptyHit();
		}    

		vec3 hitPosition = ray.origin + smallestTInInterval * ray.direction;      
		
		// To get the normal of the hit point, the point on the axis of the cylinder that normal to the hit point is needed to find. The distance of cylinder position to the
		// point on the axis can be calculted as "dot((hitPosition - cylinder.position), cylinder.direction)", since "cylinder.direction" is already a normal vector.
		// Then the result times "cylinder.direction" can get the vector of cylinder position to the "point on the axis". 
		// And add the "cylinder.position" to get the position of that point on the axis (vector + point = point).
		// Finally, get the difference of the "hitPosition" and this result, we get the normol of the hit point.
		// BTW, we need make sure the direction of the normal vector is correct and normalise it.
		vec3 normal = normalize(hitPosition - (cylinder.position + dot((hitPosition - cylinder.position), cylinder.direction) * cylinder.direction));
		
		return HitInfo
		(
		true,
		smallestTInInterval,
		hitPosition,
		normal,
		cylinder.material
		);
    }	
    return getEmptyHit();
#endif 

    return getEmptyHit();
}


HitInfo getBetterHitInfo(const HitInfo oldHitInfo, const HitInfo newHitInfo) 
{
	if(newHitInfo.hit)
  		if(newHitInfo.t < oldHitInfo.t)  // No need to test for the interval, this has to be done per-primitive
		  return newHitInfo;
  	return oldHitInfo;
}

HitInfo intersectScene(const Scene scene, const Ray ray, const float tMin, const float tMax) 
{
	HitInfo bestHitInfo;
	bestHitInfo.t = tMax;
	bestHitInfo.hit = false;
	for (int i = 0; i < cylinderCount; ++i) 
	{
		bestHitInfo = getBetterHitInfo(bestHitInfo, intersectCylinder(ray, scene.cylinders[i], tMin, tMax));
	}
	for (int i = 0; i < sphereCount; ++i) 
	{
		bestHitInfo = getBetterHitInfo(bestHitInfo, intersectSphere(ray, scene.spheres[i], tMin, tMax));
	}
	for (int i = 0; i < planeCount; ++i) 
	{
		bestHitInfo = getBetterHitInfo(bestHitInfo, intersectPlane(ray, scene.planes[i], tMin, tMax));
	}
  
	return bestHitInfo;
}

vec3 shadeFromLight
	(
	const Scene scene,
	const Ray ray,
	const HitInfo hit_info,
	const PointLight light
	)
{ 
	vec3 hitToLight = light.position - hit_info.position;
	vec3 lightDirection = normalize(hitToLight);
	vec3 viewDirection = normalize(hit_info.position - ray.origin);
	vec3 reflectedDirection = reflect(viewDirection, hit_info.normal);
	float diffuse_term = max(0.0, dot(lightDirection, hit_info.normal));
	float specular_term  = pow(max(0.0, dot(lightDirection, reflectedDirection)), hit_info.material.glossiness);

#ifdef SOLUTION_SHADOW

	// Put your shadow test here
	
	// Implement the function intersectScene() to discriminate shadow.
	// For any hit is true, the light is blocked by objects, and there is a shadow exists, 
	// then term visibility set to 0, and avoids local illumination.

	// Typical error:
	// For function intersectScene(), the fourth input tMax can not be easily set as 10000.0 or others,
	// it should be limited by length of light position to shadow hit position. 
	// Computed by function length(), so the fourth input should be length(hitToLight).
	// For some situations that object is after the light and no shadow should occur, 
	// but without this restriction, wrong shadow will be made.
	// In this scene, a little triangle-shape error shadow exists near the yellow plastic sphere if no restriction, 
	// where should not have any shadow.

	float visibility = 1.0;
	
	Ray shadowsRay;
	shadowsRay.origin = hit_info.position;
	shadowsRay.direction = lightDirection;

	HitInfo shadowsHit;
	shadowsHit = intersectScene(scene, shadowsRay, 0.001, length(hitToLight)); 

	if(shadowsHit.hit)
	{
		visibility = 0.0;
	}

#else
  	float visibility = 1.0;
#endif

	Ray mirrorRay;
	mirrorRay.origin = hit_info.position;
	mirrorRay.direction = reflect(lightDirection, hit_info.normal);
	HitInfo mirrorHitInfo = intersectScene(scene, mirrorRay, 0.001, 100000.0);
	
	return visibility * light.color * (specular_term * hit_info.material.specular + diffuse_term * hit_info.material.diffuse);
}

vec3 background(const Ray ray) 
{
	// A simple implicit sky that can be used for the background
	return vec3(0.2) + vec3(0.8, 0.6, 0.5) * max(0.0, ray.direction.y);
}

// It seems to be a WebGL issue that the third parameter needs to be inout instead of const on Tobias' machine
vec3 shade(const Scene scene, const Ray ray, inout HitInfo hitInfo) 
{
  	if(!hitInfo.hit) 
	{
		return background(ray);
  	}
  
    vec3 shading = scene.ambient * hitInfo.material.diffuse;
    for (int i = 0; i < lightCount; ++i) 
	{
		shading += shadeFromLight(scene, ray, hitInfo, scene.lights[i]); 
    }
    return shading;
}


Ray getFragCoordRay(const vec2 frag_coord) 
{
	float sensorDistance = 1.0;
	vec2 sensorMin = vec2(-1, -0.5);
	vec2 sensorMax = vec2(1, 0.5);
	vec2 pixelSize = (sensorMax- sensorMin) / vec2(800, 400);
	vec3 origin = vec3(0, 0, sensorDistance);
	vec3 direction = normalize(vec3(sensorMin + pixelSize * frag_coord, -sensorDistance));  
  
  	return Ray(origin, direction);
}

float fresnel(const vec3 viewDirection, const vec3 normal) 
{
#ifdef SOLUTION_FRESNEL
	// Put your code to compute the Fresnel effect here
	
	// When a light have refraction from one medium to the other medium, the reflection may also occur, 
	// the Fresnel factor presents the weighting relationship between refraction and reflection.
	// For Fresnel effect, the Schlick's approximation can be applied to approximate the contribution of 
	// the Fresnel factor in a specular reflection between two surface.

	// In Schlick's model, the approximation coefficient of reflection can be presentd by R(theta):
	// R(theta) = R0 + (1 + R0) * (1 - cos(theta))^5
	// R0 = ((n1 - n2)/(n1 + n2))^2
	
	// Where, theta is angle between the vector of ray direction and normal direction,
	// and n1 and n2 are indices of refraction of two mediums,
	// n1 is the IOR of next material and n2 is the IOR of previous material. 
	// In this coursework, the IOR of air is 1.0 and the IOR of glass is 1.05 for best result.
	// Then R0 can be pre-set as 0.00059, ((1.05 - 1)/(1.05 + 1))^2 = 0.00059
	// For refraction, the approximation coefficient can be presented by: (1 - R(theta)).
	// Since we assume W_reflect = 1 - W_refract.
	// Additional, for best fitting the result from coursework instruction, the power coefficient set by 1.6

	float R0 = 0.00059;
	float cosTheta = dot(-viewDirection, normal);
	float rTheta = R0 + (1.0 - R0) * pow((1.0 - cosTheta), 1.6);
	return rTheta;
#else
	return 1.0;
#endif
}

vec3 colorForFragment(const Scene scene, const vec2 fragCoord) 
{
      
  	Ray initialRay = getFragCoordRay(fragCoord);  
  	HitInfo initialHitInfo = intersectScene(scene, initialRay, 0.001, 10000.0);  
  	vec3 result = shade(scene, initialRay, initialHitInfo);
	
  	Ray currentRay;
  	HitInfo currentHitInfo;
  	
  	// Compute the reflection
  	currentRay = initialRay;
  	currentHitInfo = initialHitInfo;
  	
  	// The initial strength of the reflection
  	float reflectionWeight = 1.0;
  	
  	const int maxReflectionStepCount = 2;
  	for (int i = 0; i < maxReflectionStepCount; i++) 
	{
		if (!currentHitInfo.hit) break;
      
#ifdef SOLUTION_REFLECTION_REFRACTION
		// Put your reflection weighting code here

		// Since each material has its own reflection ability.  
		reflectionWeight *= currentHitInfo.material.reflection;    
#endif
      
#ifdef SOLUTION_FRESNEL
		// Add Fresnel contribution

		// The reflection strength was influenced by the Fresnel contribution.
		reflectionWeight = reflectionWeight * fresnel(currentRay.direction, currentHitInfo.normal);
#else
		reflectionWeight *= 0.5;
#endif
      
		Ray nextRay;
#ifdef SOLUTION_REFLECTION_REFRACTION
		// Put your code to compute the reflection ray here

		// For each reflection, the incident angle is equal to the reflection angle,
		// the incident angle computed by incident ray to the normal,
		// and the reflection angle computed by the reflection ray to the normal.
		// And the next ray continues from the hit position.
		nextRay.direction = reflect(currentRay.direction, currentHitInfo.normal);
		nextRay.origin = currentHitInfo.position;
#endif
		currentRay = nextRay; 
		currentHitInfo = intersectScene(scene, currentRay, 0.001, 10000.0);             
		result += reflectionWeight * shade(scene, currentRay, currentHitInfo);
	}
  
	// Compute the refraction
	currentRay = initialRay;  
	currentHitInfo = initialHitInfo;
   
  	// The initial medium is air
  	float currentIOR = 1.0;

  	// The initial strength of the refraction.
  	float refractionWeight = 1.0;
  
  	const int maxRefractionStepCount = 2;
  	for(int i = 0; i < maxRefractionStepCount; i++) 
	{
    
#ifdef SOLUTION_REFLECTION_REFRACTION
		// Put your refraction weighting code here

		// Materials have their own refraction ability.
		// And for opaque materials, these terms will be set as 0,
		// then break the refraction calculation loop.
		refractionWeight *= currentHitInfo.material.refraction;
		if (currentHitInfo.material.refraction == 0.0) break;
#else
		refractionWeight *= 0.5;      
#endif

#ifdef SOLUTION_FRESNEL
		// Add Fresnel contribution

		// The refraction strength was influenced by the Fresnel contribution.
		// The sum of coefficients of reflection and refraction equal to 1.
		refractionWeight = refractionWeight * (1.0 - fresnel(currentRay.direction, currentHitInfo.normal));
#endif      

		Ray nextRay;
#ifdef SOLUTION_REFLECTION_REFRACTION      
		// Put your code to compute the refraction ray and track the IOR
		
		// Since in different medium light travels with varies velocities,
		// refraction occurs between two materials.
		// By Snell's Law, 
		// n_1 * Sin(theta_1) = n_2 * Sin(theta_2)
		// Where, n_i is the index of refraction of materials,
		// and theta_i is the angle between the light and the normal.

		// Implement build-in function refract() to calculate the refraction ray.
		// Input ray direction, normal direction and IOR ratio to get the refraction ray.
		// Since this scene is not so complex, to reduce the computing complexity,
		// we can force the nextIOR of second iteration is 1,
		// which means when ray comes out from glass sphere to the air, the IOR of air is 1.
		
		float nextIOR;
		nextIOR = currentHitInfo.material.IOR;
		if(i==1){nextIOR = 1.0;}
		
		nextRay.direction = refract(currentRay.direction, currentHitInfo.normal, currentIOR/nextIOR);
		nextRay.origin = currentHitInfo.position;
		currentIOR = nextIOR;

#endif
		currentRay = nextRay;
		currentHitInfo = intersectScene(scene, currentRay, 0.001, 10000.0);   
		result += refractionWeight * shade(scene, currentRay, currentHitInfo);
      
		if (!currentHitInfo.hit) break;
	}
	return result;
}

Material getDefaultMaterial() 
{
#ifdef SOLUTION_MATERIAL
	// Update the default material call to match the new parameters of Material
	return Material(vec3(0.3), vec3(0), 1.0, 0.0, 0.0, 1.0);
#else
	return Material(vec3(0.3), vec3(0), 1.0);
#endif
}
	// Explanation of terms in Material struct:

	// Diffuse: the diffuse reflection generally occurs on a roughness surface,
	// and the reflection lights will randomly reflect to many directions.

	// Specular: the specular reflection occurs on a smooth surface,
	// each reflection light has a same angle to normal with the incident light.

	// Glossiness: Gloss presents how incident lights reflect on a specular surface,
	// a higher glossiness leads a smaller highlight area on a surface, and vice versa.

	// reflection: reflection ability of this material.

	// refraction: refraction ability of this material, set it as 0 for opaque material.

	// IOR: Index of refraction.


Material getPaperMaterial() 
{
#ifdef SOLUTION_MATERIAL
	// Definition of a paper material
	return Material(
		vec3(0.95), 				// diffuse, paper has a roughness surface, and for a white paper,
									// each channel should be very high.
		vec3(0.1), 					// specular, and so the specular should be low.
		2.0, 						// glossiness, glossiness is low.
		0.1,						// reflection, reflection ability of paper is low.
		0.0,						// refraction, paper is opaque.
		0.0							// IOR
		);
#else
    return getDefaultMaterial();
#endif
}

Material getPlasticMaterial() 
{
#ifdef SOLUTION_MATERIAL
	// Definition of a plastic material
	return Material(
		vec3(0.5, 0.2, 0.0),		// diffuse, for a yellow or orange plastic, the channel R should be higher,
									// and with a bit contribution in channel G, to make the material looks like orange.
		vec3(0.7),					// specular, plastic has a smooth surface, specular should be a bit higher.
		8.0,						// glossiness, glossiness of plastic is higher than paper but much lower than mirror.
		0.6,						// reflection, good reflection ability, but lower than mirror.
		0.0,						// refraction, plastic is opaque.
		0.0							// IOR
		);
#else
    return getDefaultMaterial();
#endif
}

Material getGlassMaterial() 
{
#ifdef SOLUTION_MATERIAL
	// Replace by your definition of a glass material
	return Material(
		vec3(0.0),					// diffuse, As for a transparent material, diffuse is zeros for no colour.
		vec3(0.0), 					// specular, specular is also zeros.
		1.0,						// glossiness, Since diffuse and specular are zeros, glossiness just set as non-zero.
		1.0,						// reflection, Glass has a good reflection ability.
		1.0,						// refraction, Glass as a good refraction ability.
		1.05						// IOR, IOR of glass, 1.05 best fitting the sample picture.
		);
#else
    return getDefaultMaterial();
#endif
}

Material getSteelMirrorMaterial() 
{
#ifdef SOLUTION_MATERIAL
	// Definition of a steel mirror material
	return Material(
		vec3(0.1),					// diffuse, For a metallic ground, it should looks gray,
									// and so all three channels should be low.
		vec3(0.9),					// specular, A steel mirror should has a high specular value.
		100.0,						// glossiness, a polished matel has a high glossiness.
		1.0,						// reflection, A mirror has a great reflection ability.
		0.0,						// refraction, steel is opaque.
		0.0							// IOR
		);
#else
    return getDefaultMaterial();
#endif
}

vec3 tonemap(const vec3 radiance) {
	const float monitorGamma = 2.0;
	return pow(radiance, vec3(1.0 / monitorGamma));
}

void main()
{
    // Setup scene
	Scene scene;
	scene.ambient = vec3(0.12, 0.15, 0.2);
  
    // Lights
	scene.lights[0].position = vec3(5, 15, -5);
	scene.lights[0].color    = 0.5 * vec3(0.9, 0.5, 0.1);
    
	scene.lights[1].position = vec3(-15, 5, 2);
	scene.lights[1].color    = 0.5 * vec3(0.1, 0.3, 1.0);
  
    // Primitives
	scene.spheres[0].position		= vec3(10, -5, -16);
	scene.spheres[0].radius			= 6.0;
	scene.spheres[0].material		= getPaperMaterial();
    
	scene.spheres[1].position		= vec3(-7, -1, -13);
	scene.spheres[1].radius			= 4.0;
	scene.spheres[1].material		= getPlasticMaterial();
  
	scene.spheres[2].position		= vec3(0, 0.5, -5);
	scene.spheres[2].radius			= 2.0;
	scene.spheres[2].material		= getGlassMaterial();

	scene.planes[0].normal			= vec3(0, 1, 0);
	scene.planes[0].d				= 4.5;
	scene.planes[0].material		= getSteelMirrorMaterial();

	scene.cylinders[0].position		= vec3(-1, 1, -18);
	scene.cylinders[0].direction	= normalize(vec3(-1, 2, -1));
	scene.cylinders[0].radius		= 1.5;
	scene.cylinders[0].material		= getPaperMaterial();
  
	scene.cylinders[1].position		= vec3(4, 1, -5);
	scene.cylinders[1].direction	= normalize(vec3(1, 4, 1));
	scene.cylinders[1].radius		= 0.4;
	scene.cylinders[1].material		= getPlasticMaterial();

	// compute color for fragment
	gl_FragColor.rgb = tonemap(colorForFragment(scene, gl_FragCoord.xy));
	gl_FragColor.a = 1.0;

}