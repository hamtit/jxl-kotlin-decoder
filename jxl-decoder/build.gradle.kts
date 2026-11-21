plugins {
    kotlin("jvm")
    `java-library`
    `maven-publish`
}

kotlin {
    // Android 12 (API 31) 以降の ART でも動くように Java 17 のバイトコードにする
    compilerOptions {
        jvmTarget.set(org.jetbrains.kotlin.gradle.dsl.JvmTarget.JVM_17)
        // CI の lint では -Pwerror を付けて警告をエラーにする
        allWarningsAsErrors.set(providers.gradleProperty("werror").isPresent)
    }
    explicitApi()
}

java {
    sourceCompatibility = JavaVersion.VERSION_17
    targetCompatibility = JavaVersion.VERSION_17
    withSourcesJar()
}

dependencies {
    testImplementation(kotlin("test"))
}

// ライセンスを jar (クラスとソースの両方) に入れる
tasks.withType<Jar>().configureEach {
    from(rootProject.files("LICENSE", "LICENSE-jxl-oxide")) { into("META-INF") }
}

tasks.test {
    useJUnitPlatform()
    maxHeapSize = "2g"
    systemProperty("jxl.testdata", rootProject.file("testdata").absolutePath)
}

publishing {
    publications {
        create<MavenPublication>("maven") {
            artifactId = "jxl-decoder"
            from(components["java"])
        }
    }
}

